export TransferFunction, sampleTF, setTF, combine

"""
$(TYPEDEF)
$(TYPEDFIELDS)


    TransferFunction(freq::Vector, data::Array, [phasedata]; [inductionFactor], [units])

Create a `TransferFunction` from a data array `data` at frequencies `freq`.
`freq` is given in Hz or should have a Unitful frequency unit attached
`data` should have the shape [frequencies, channels]. `data` can be either complex or real, if `data` is real a second array `phasedata` can be passed representing the phase in radians. Both the amplitude and the phase can have Unitful units.

# Optional Keyword-Arguments:
- `inductionFactor::Vector{<:Real}`: induction factor for each channel
- `units::Vector`: units for each channel, can be either Unitful.FreeUnits or a string that can be parsed as a Unitful unit. Instead of using this keyword, `data` can also have Unitful units attached, then the units keyword-argument is ignored.

"""
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

function TransferFunction(freq::Vector{<:Real}, ampdata::Array{<:Union{Real,Unitful.AbstractQuantity{<:Real}},N}, phasedata::Array{<:Real,N}; kwargs...) where N
  if size(ampdata) != size(phasedata); error("The size of ampdata and phasedata must match!") end
  data = ampdata.*exp.(im.*phasedata)
  return TransferFunction(freq, data; kwargs...)
end

TransferFunction(freq::Vector{<:Unitful.Frequency}, args...; kwargs...) = TransferFunction(ustrip.(u"Hz", freq), args...; kwargs...)
TransferFunction(freq::Vector{<:Real}, ampdata::Array{<:Union{Real,Unitful.AbstractQuantity{<:Real}},N}, phasedata::Array{<:Unitful.DimensionlessQuantity{<:Real},N}; kwargs...) where N = TransferFunction(freq, ampdata, ustrip.(u"rad", phasedata); kwargs...)
TransferFunction(freq::Vector{<:Real}, ampdata::Array{<:Union{Real,Unitful.AbstractQuantity{<:Real}},N}; kwargs...) where N = TransferFunction(freq, ampdata, zeros(size(ampdata)); kwargs...)

function TransferFunction(freq::Vector{<:Real}, datain::Array{<:Unitful.AbstractQuantity{<:Complex},N}; units=nothing, kwargs...) where N
  if !isnothing(units)
    @warn "You passed explicit units combined with a Unitful data array. The explicit units will be ignored!"
  end
  units = Unitful.FreeUnits[]
  for ch in eachcol(datain)
    ch_units = Unitful.unit.(ch)
    if all(y->y==ch_units[1], ch_units)
      push!(units, ch_units[1])
    else
      error("One TransferFunction channel must have identical units!")
    end
  end
  return TransferFunction(freq, ustrip.(datain); units=units, kwargs...)
end

"""
$(TYPEDSIGNATURES)

Create a `TransferFunction` from a data file at `filename`.

The file can be either a h5-File created with this package or a file that is written by a VNA. Keyword arguments will be passed to `load_tf_fromVNA` 
"""
function TransferFunction(filename::String; kargs...)
    filenamebase, ext = splitext(filename)
    if ext == ".h5"
      tf = load_tf(filename; kargs...)
    else #if ext == "s1p" || ext == "s2p"
      tf = load_tf_fromVNA(filename; kargs...)
    end
    return tf
end

"""
$(TYPEDSIGNATURES)

Create a `TransferFunction` from the tf data saved in a MPIFile (see `rxTransferFunction`)
"""
function TransferFunction(file::MPIFile)
  tf_file = rxTransferFunction(file)
  if isnothing(tf_file)
    return nothing
  end
  inductionFactor = rxInductionFactor(file)
  f = rxFrequencies(file)
  
  if isnothing(inductionFactor)
    return TransferFunction(f, abs.(tf_file), angle.(tf_file))
  else
    return TransferFunction(f, abs.(tf_file), angle.(tf_file), inductionFactor=inductionFactor)
  end
end

"""
    tf[i,j]

Directly access the underlying data of a `TransferFunction`
"""
function getindex(tf::TransferFunction, args...)
  try 
    return getindex(tf.data, args...)  
  catch e
    @warn "The indexing using square brackets on TransferFunction objects now always operates on the integer indizes of the underlying transfer function data. To use frequency interpolation, use tf(freq, channel) instead of tf[[freq],channel]."
    rethrow(e)
  end
end

"""
    tf(f, chan::Integer)

Interpolated access to a `TransferFunction` at frequencies `f` and single channel `chan`
"""
function (tf::TransferFunction)(f, chan::Integer=1)
  if chan>length(tf.interpolator); error("The TransferFunction only has $(length(tf.interpolator)) channel(s), unable to access channel $(chan)") end
  return tf.interpolator[chan](f) .* tf.units[chan]
end

"""
    tf(f, chan::AbstractVector{<:Integer})

Interpolated access to a `TransferFunction` at frequencies `f` and channels `chan`
"""
(tf::TransferFunction)(f, chan::AbstractVector{<:Integer}) = hcat([tf(f,c) for c in chan]...)

(tf::TransferFunction)(f, ::Colon) = tf(f, axes(tf.data,2))

"""
    load_tf(filename::String; channels = nothing, kwargs...)

Load a `TransferFunction` from a h5-file at `filename`. If `channels` is set, only the specified channels are loaded.
"""
function load_tf(filename::String; channels = nothing)
  if channels isa Integer
    channels = [channels]
  end

  return h5open(filename, "r") do file
    tf = read(file, "/transferFunction")
    freq = read(file,"/frequencies")
    inductionFactor = read(file,"/inductionFactor")

    if !isnothing(channels)
      tf = tf[:,channels]
      inductionFactor = inductionFactor[channels]
    end

    if haskey(file, "/unit")
      unit = uparse.(read(file, "/unit"))

      if !isnothing(channels)
        unit = unit[channels]
      end

      return TransferFunction(freq, tf, inductionFactor=inductionFactor, units=unit)
    end

    return TransferFunction(freq, tf, inductionFactor=inductionFactor)
  end
end

"""
$(TYPEDSIGNATURES)
Combine two `TransferFunctions` along their channel dimension. If interpolate=false, will only work if the frequency samples are identical.
"""
function combine(tf1::TransferFunction, tf2::TransferFunction; interpolate=false)
  if !interpolate
    if tf1.freq != tf2.freq; error("The frequency axes of the transfer functions do not match. Can not combine!") end
    freq = tf1.freq
    data = cat(tf1.data,tf2.data, dims=2)
  else
    freq = unique(sort(cat(tf1.freq, tf2.freq, dims=1)))
    data = cat(ustrip.(tf1(freq, :)), ustrip.(tf2(freq, :)), dims=2)
  end
    inductionFactor = cat(tf1.inductionFactor, tf2.inductionFactor, dims=1)
    units = cat(tf1.units, tf2.units, dims=1)
    return TransferFunction(freq, data, inductionFactor=inductionFactor, units=units)
end

combine(tf::TransferFunction, ::Nothing; interpolate=false) = combine(tf, TransferFunction(tf.freq, zeros(length(tf.freq))))
combine(::Nothing, tf::TransferFunction; interpolate=false) = combine(TransferFunction(tf.freq, zeros(length(tf.freq))), tf) 

"""
$(TYPEDSIGNATURES)
Save `tf` as a h5 file to `filename`
"""
function save(filename::String, tf::TransferFunction)
  h5write(filename, "/transferFunction", tf.data)
  h5write(filename, "/frequencies", tf.freq)
  h5write(filename, "/inductionFactor", tf.inductionFactor)
  h5write(filename, "/unit", string.(tf.units))
  return nothing
end

"""
$(TYPEDSIGNATURES)
Read the data recorded with the VNA from file `filename` into three Julia arrays `freq`, `ampdata`, `phasedata`
"""
function readVNAdata(filename::String)
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
  return freq, apdata, aϕdata
end

"""
$(TYPEDSIGNATURES)
Load data receive calibration from file recorded with the VNA. Data will be processed by `processRxTransferFunction`, see there for keyword arguments
"""
function load_tf_fromVNA(filename::String; kwargs...)
  freq, ampdata, phasedata = readVNAdata(filename)
  compdata = ampdata.*exp.(im.*phasedata)
  if any([:R,:N,:A] .∈ [keys(kwargs)])
    @warn "In v0.16.1 and below load_tf_fromVNA mistakenly ignored the keyword parameters R, N and A. The current version includes these parametes if set, resulting in the correct scaling in magnetic moment domain."
  end
  return processRxTransferFunction(freq, compdata; kwargs...)    
end

"""
$(TYPEDSIGNATURES)
Process the data from a receive calibration measurement using a calibration coil into a TransferFunction

Keyword parameters:
- `frequencyWeighting`: if true corrects for the frequency term in the TF, which results in a TF that does not integrate on application, but instead shows derivative of magnetic moment
- `R`: value of the resistance of the calibration coil in Ω
- `N`: number of turns of the calibration coil
- `A`, `r`, `d`: Area in m² or radius/diameter in m of the calibration coil, define only one
"""
function processRxTransferFunction(freq, compdata; frequencyWeighting::Bool=false,
  R::Union{Real,Nothing} = nothing, # Ω
  N::Union{Real,Nothing} = nothing, # Turns
  A::Union{Real,Nothing} = nothing, # m²
  r::Union{Real,Nothing} = nothing, # m
  d::Union{Real,Nothing} = nothing) # m

  if sum(.!isnothing.([A,r,d])) > 1
    error("You can only define one of the keyword parameters A, r or d defining the geometry of the calibration coil")
  elseif !isnothing(d)
    A = pi * (d/2)^2
  elseif !isnothing(r)
    A = pi * r^2
  end

  if all(isnothing.([R,N,A]))
    unit = u"V/V"
    @warn "No calibration coil parameters set in `processRxTransferFunction`. Make sure this is intended"
  elseif any(isnothing.([R,N,A]))
    error("To process the rxTransferFunction all three parameters describing the calibration coil (R, N and A) need to be defined.")
  else
    unit = u"V/A*m^2"
    compdata .*= R/(N*A)
  end
  
  if frequencyWeighting
    compdata ./= (im*freq.*2*pi) # As TF is defined as u_ADC = u_coil *TF the derivative from magnetic moment is applied as the division of the TF by jw
  end

  return TransferFunction(freq, compdata, units=[unit])
end

"""
$(TYPEDSIGNATURES)
Sample the `tf` at the frequencies defined by the measurement `f`, optionally with `numPeriodGrouping`
"""
function sampleTF(tf::TransferFunction, f::MPIFile; numPeriodGrouping=1)
  freq = rxFrequencies(f, numPeriodGrouping)
  if measIsFrequencySelection(f)
    freq = freq[measFrequencySelection(f)]
  end
  numChan = rxNumChannels(f)
  return tf(freq,1:numChan)
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

  if isnothing(l)
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
  if !isnothing(g)
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

setTF(f::MDFFile, filenameTF::AbstractString) = setTF(f, TransferFunction(filenameTF), filenameTF)

function setTF(f::MDFFile, tmf::TransferFunction, filenameTF::Union{AbstractString,Nothing}=nothing)
  tf = sampleTF(tmf, f)

  # We need to close the HDF5 file handle before we can write to it
  close(f.file)

  h5open(filepath(f), "r+") do file
	  if haskey(file, "/acquisition/receiver/transferFunctionFileName")
      delete_object(file, "/acquisition/receiver/transferFunctionFileName")
    end
    if !isnothing(filenameTF)
      write(file, "/acquisition/receiver/transferFunctionFileName", filenameTF)
    end
    if haskey(file, "/acquisition/receiver/transferFunction")
      delete_object(file, "/acquisition/receiver/transferFunction")
    end
    write(file, "/acquisition/receiver/transferFunction", tf)
    if haskey(file, "/acquisition/receiver/inductionFactor")
      delete_object(file, "/acquisition/receiver/inductionFactor")
    end
    write(file, "/acquisition/receiver/inductionFactor", tmf.inductionFactor)
  end
  # reopen filehandler so f stays usable
  f.file = h5open(filepath(f), "r")
  return
end

function setTF(f::MDFv2InMemory, tf::TransferFunction)
  rxTransferFunction(f, sampleTF(tf, f))
  rxInductionFactor(f, tf.inductionFactor)
  return
end