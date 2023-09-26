export TransferFunction, sampleTF, setTF

mutable struct TransferFunction
  freq::Vector{Float64}
  data::Matrix{ComplexF64}
  inductionFactor::Vector{Float64}

  function TransferFunction(freq_, datain::Array{T}, inductionFactor=ones(size(datain, 2))) where {T<:Complex}
    freq = freq_[1]:(freq_[2]-freq_[1]):freq_[end]
    data = reshape(deepcopy(datain), size(datain,1), size(datain, 2))
    return new(freq_, data, inductionFactor)
  end
end

function TransferFunction(freq_, ampdata, phasedata, args...)
  data = ampdata.*exp.(im.*phasedata)
  return TransferFunction(freq_, data, args...)
end

function TransferFunction(filename::String; kargs...)
    filenamebase, ext = splitext(filename)
    if  ext == ".h5"
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

function getindex(tmf::TransferFunction, x::UnitRange, chan::Integer=1)
  a = tmf.data[x,chan]
  return a
end

#function getindex(tmf::TransferFunction, x::Real, chan::Integer=1)
#  a = tmf.interp[chan](x)
#  return a
#end

function getindex(tmf::TransferFunction, X::Union{Vector,AbstractRange}, chan::Integer=1)
  I = extrapolate(interpolate((tmf.freq,), tmf.data[:,chan], Gridded(Linear())), Interpolations.Flat())

  return [I(x) for x in X]
end

function getindex(tmf::TransferFunction, X::Union{Vector,AbstractRange}, chan::AbstractRange)
  out = zeros(ComplexF64, length(X), length(chan))
  for d=1:length(chan)
    out[:,d] = tmf[X,d]
  end
  return out
end

function load_tf(filename::String)
  tf = h5read(filename,"/transferFunction")
  freq = h5read(filename,"/frequencies")
  inductionFactor = h5read(filename,"/inductionFactor")
  return TransferFunction(freq,tf,inductionFactor)
end

function combine(tf1, tf2)
  freq = tf1.freq
  data = cat(tf1.data,tf2.data, dims=2)
  inductionFactor = cat(tf1.inductionFactor, tf2.inductionFactor, dims=1)
  return TransferFunction(freq, data, inductionFactor)
end

function save(filename::String, tf::TransferFunction)
  h5write(filename, "/transferFunction", tf.data)
  h5write(filename, "/frequencies", tf.freq)
  h5write(filename, "/inductionFactor", tf.inductionFactor)
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
  line_trafo = nothing
  for lineIdx in eachindex(lines)
    line = lines[lineIdx]
    if startswith(line, "!") # This is a comment and thus we ignore this line
      continue
    elseif startswith(line, "#")
      header = split(line, " ", keepempty=false)
      frequencyUnit = header[2]
      parameter = header[3] # Assumed to be S for "Scattering parameters"
      mode = header[4]
      impedance = header[6] # header[4] is always "R"

      frequencyFactor = 1
      if lowercase(frequencyUnit) == "hz"
        frequencyFactor = 1
      elseif lowercase(frequencyUnit) == "khz"
        frequencyFactor = 1e3
      elseif lowercase(frequencyUnit) == "mhz"
        frequencyFactor = 1e6
      elseif lowercase(frequencyUnit) == "ghz"
        frequencyFactor = 1e9
      else
        error("Frequency unit `$frequencyUnit` not applicable.")
      end

      line_trafo = (line) -> begin
        line_parts = split(line, " ", keepempty=false)
        f = parse(Float64, line_parts[1])*frequencyFactor

        if uppercase(mode) == "DB" # dB-angle (dB = 20*log10|magnitude|)
          ap = 10^(parse(Float64, line_parts[2])/20)
          aphi = parse(Float64, line_parts[3])*π/180
        elseif uppercase(mode) == "MA" # magnitude-angle
          ap = parse(Float64, line_parts[2])
          aphi = parse(Float64, line_parts[3])
        elseif uppercase(mode) == "RI" # real-imaginary
          tf_complex=parse(Float64, line_parts[2]) + im*parse(Float64, line_parts[3])
          ap = abs(tf_complex)
          aphi = angle(tf_complex)
        else
          error("Mode `$mode` not applicable.")
        end

        return ap, aphi, f
      end
    else
      if !isnothing(line_trafo)
        ap, aphi, f = line_trafo(line)
        push!(apdata, ap)
        push!(aϕdata, aphi)
        push!(freq, f)
      else
        error("No header has been found until first data lines were reached.")
      end
    end
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
  numFreq = length(freq)
  return tmf[freq,1:numChan]
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
