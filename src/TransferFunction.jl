export TransferFunction, sampleTF, setTF

mutable struct TransferFunction
  freq::Vector{Float64}
  data::Matrix{ComplexF64}
  interpolator::Vector{AbstractInterpolation}
  inductionFactor::Vector{Float64}
  units::Vector{Unitful.FreeUnits}

  function TransferFunction(freq_::Vector{<:Real}, datain::Array{<:Complex}; inductionFactor::Vector{<:Real}=ones(size(datain, 2)), units::Vector=Unitful.FreeUnits[Unitful.NoUnits for i in 1:size(datain, 2)])

    parsed_units = Unitful.FreeUnits[]
    for tmp in units
      if isa(tmp, String); tmp = uparse(tmp) end # get correct unit from unit strings
      if !isa(tmp, Unitful.FreeUnits); tmp = Unitful.unit(tmp) end # get correct unit for numbers and quantities, e.g. unit(1) or unit(1u"V/T"), includes the case that the string parsing returned one of these types
      push!(parsed_units, tmp)
    end
    
    if length(freq_) != size(datain,1); error("Length of frequency axis ($(length(freq_))) does not match the number of samples in the data ($(size(datain,1)))!") end
    if size(datain, 2) != length(inductionFactor); error("Number of channels in data ($(size(datain, 2))) does not match the number of channels in the given inductionFactor ($(size(inductionFactor,1)))!") end
    if size(datain, 2) != length(units); error("Number of channels in data ($(size(datain, 2))) does not match the number of channels in the given units ($(size(units,1)))!") end

    data = reshape(deepcopy(datain), size(datain,1), size(datain, 2))
    interpolator = [extrapolate(interpolate((freq_,), channel, Gridded(Linear())), Interpolations.Flat()) for channel in eachcol(data)]
    return new(freq_, data, interpolator, inductionFactor, parsed_units)
  end
end

Base.show(io::IO, ::MIME"text/plain", tf::TransferFunction) = print(io, "MPIFiles.TransferFunction: \n\t$(size(tf.data,2)) channel(s), units of $(string.(tf.units))\n\t$(size(tf.data,1)) frequency samples from $(tf.freq[1]) Hz to $(tf.freq[end]) Hz")

function TransferFunction(freq_::Vector{<:Real}, ampdata::Array{<:Real,N}, phasedata::Array{<:Real,N}; kwargs...) where N
  if size(ampdata) != size(phasedata); error("The size of ampdata and phasedata must match!") end
  data = ampdata.*exp.(im.*phasedata)
  return TransferFunction(freq_, data; kwargs...)
end

function TransferFunction(filename::String; kargs...)
    filenamebase, ext = splitext(filename)
    if ext == ".h5"
      tf = load_tf(filename)
    else #if ext == "s1p" || ext == "s2p"
      tf = load_tf_fromVNA(filename; kargs...)
    end
    return tf
end

function TransferFunction(file::MPIFile)
  tf_file = rxTransferFunction(file)
  inductionFactor = rxInductionFactor(file)
  f = collect(rfftfreq(rxNumSamplingPoints(file), rxBandwidth(file)*2))
  return TransferFunction(f, abs.(tf_file), angle.(tf_file), inductionFactor)
end

function getindex(tmf::TransferFunction, args...)
  try 
    return getindex(tmf.data, args...)  
  catch e
    @warn "The indexing using square brackets on TransferFunction objects now always operates on the integer indizes of the underlying transfer function data. To use frequency interpolation, use tf(freq, channel) instead of tf[[freq],channel]."
    rethrow(e)
  end
end

function (tmf::TransferFunction)(x, chan::Integer=1)
  if chan>length(tmf.interpolator); error("The TransferFunction only has $(length(tmf.interpolator)) channel(s), unable to access channel $(chan)") end
  return tmf.interpolator[chan](x) .* tmf.units[chan]
end

(tmf::TransferFunction)(x, chan::AbstractArray) = hcat([tmf(x,c) for c in chan]...)

function load_tf(filename::String)
  tf = h5read(filename,"/transferFunction")
  freq = h5read(filename,"/frequencies")
  inductionFactor = h5read(filename,"/inductionFactor")
  unit = []
  try 
    unit = h5read(filename, "/unit")
  catch # if h5read fails, it should mean that there is no units, maybe do something better here
    return TransferFunction(freq, tf, inductionFactor=inductionFactor)
  end
  return TransferFunction(freq, tf, inductionFactor=inductionFactor, units=uparse.(unit))
end

function combine(tf1::TransferFunction, tf2::TransferFunction)
  if tf1.freq != tf2.freq; error("The frequency axes of the transfer functions do not match. Can not combine!") end
  freq = tf1.freq
  data = cat(tf1.data,tf2.data, dims=2)
  inductionFactor = cat(tf1.inductionFactor, tf2.inductionFactor, dims=1)
  units = cat(tf1.units, tf2.units, dims=1)
  return TransferFunction(freq, data, inductionFactor=inductionFactor, units=units)
end

function save(filename::String, tf::TransferFunction)
  h5write(filename, "/transferFunction", tf.data)
  h5write(filename, "/frequencies", tf.freq)
  h5write(filename, "/inductionFactor", tf.inductionFactor)
  h5write(filename, "/unit", string.(tf.units))
  return nothing
end

function load_tf_fromVNA(filename::String;
    frequencyWeighting=false,
    R = 50.0, #Ω
    N = 10, #5# Turns
    A = 7.4894*10.0^-4) #m^2 #1.3e-3^2*pi;)

  file = open(filename)
  lines = readlines(file)
  apdata = Float64[]
  aϕdata = Float64[]
  freq = Float64[]
  for i=4:length(lines)
      tmp = split(strip(lines[i])," ")
      tmp=tmp[tmp.!=""]
      f = parse(Float64,strip(tmp[1]))
      if occursin(lines[3],"# kHz S MA R 50")
          ap = parse(Float64,strip(tmp[2]))
          aphi = parse(Float64,strip(tmp[3]))
          f=f*1000
       elseif occursin(lines[3],"# Hz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f
       elseif occursin(lines[3],"# kHz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f*1000
        elseif occursin(lines[3],"# MHz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f*1000000
      elseif occursin(lines[3],"# kHz S RI R 50")
          tf_complex=parse(Float64,strip(tmp[2]))+im*parse(Float64,strip(tmp[3]))
          ap=abs.(tf_complex);
          aphi=angle(tf_complex)
          f=f*1000
      else
	      error("Wrong data Format! Please export in kHz domain S21 parameter with either Magnitude/Phase, DB/Phase or Real/Imaginary!")
      end
      push!(apdata, ap)
      push!(aϕdata, aphi)
      push!(freq, f)
  end
  close(file)
  if frequencyWeighting
  	apdata ./= (freq.*2*pi) # As TF is defined as u_ADC = u_coil *TF the derivative from magnetic moment is applied as the division of the TF by w
  end
  return TransferFunction(freq, apdata, aϕdata)
end


function sampleTF(tmf::TransferFunction, f::MPIFile)
  freq = rxFrequencies(f)
  numChan = rxNumChannels(f)
  return tmf(freq,1:numChan)
end


function _writeBrukerParamFile(b::BrukerFile, A, filename)
  mv(joinpath(filepath(b),filename), joinpath(filepath(b),filename*"_bac"), force=true )

  open(joinpath(filepath(b),filename),"w") do fd
    for u in A
      println(fd,u)
    end
  end
  return
end

function setTF(b::BrukerFile, filenameTF::AbstractString)
  A = readlines(joinpath(filepath(b),"acqp"))
  l = findfirst( s->occursin("ACQ_comment", s), A)

  if l == nothing
   ll = findfirst( s->occursin("ACQ_experiment_mode", s), A)
   insert!(A, ll, "##\$ACQ_comment=( 2048 )")
   insert!(A, ll+1, "<$(filenameTF)>")
  else
    A[l+1] = "<$(filenameTF)>"
  end

  _writeBrukerParamFile(b, A, "acqp")

  # For some reason ACQ_comment is also in the method file -> also change that
  B = readlines(joinpath(filepath(b),"method"))
  g = findfirst( s->occursin("ACQ_comment", s), B)
  if g != nothing
    B[g+1] = "<$(filenameTF)>"
    _writeBrukerParamFile(b, B, "method")
  end

  if isfile(joinpath(filepath(b),"mdf"))
    MPIFile(filepath(b)) do f
      setTF(f, filenameTF)
    end
  end
  return
end


function setTF(f::MDFFile, filenameTF::AbstractString)
  tmf = TransferFunction(filenameTF)
  tf = sampleTF(tmf, f)

  # We need to close the HDF5 file handle before we can write to it
  close(f.file)

  h5open(filepath(f), "r+") do file
	  if haskey(file, "/acquisition/receiver/transferFunctionFileName")
      delete_object(file, "/acquisition/receiver/transferFunctionFileName")
    end
    write(file, "/acquisition/receiver/transferFunctionFileName", filenameTF)
    if haskey(file, "/acquisition/receiver/transferFunction")
      delete_object(file, "/acquisition/receiver/transferFunction")
    end
    write(file, "/acquisition/receiver/transferFunction", tf)
    if haskey(file, "/acquisition/receiver/inductionFactor")
      delete_object(file, "/acquisition/receiver/inductionFactor")
    end
    write(file, "/acquisition/receiver/inductionFactor", tmf.inductionFactor)
  end
  return
end
