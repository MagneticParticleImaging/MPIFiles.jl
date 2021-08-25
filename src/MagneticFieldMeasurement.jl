export MagneticFieldMeasurement, saveMagneticFieldAsHDF5, loadMagneticFieldMeasurement, loadMagneticField

Base.@kwdef mutable struct MagneticFieldMeasurement
  "Description of the dataset."
  description::Union{String, Nothing} = nothing
  "Positions of the measured field values."
  positions::Union{Positions, Missing} = missing
  "Field values at the specific positions."
  fields::Union{Array{typeof(1.0u"T"), 2}, Missing} = missing
  "Error of the measured field values."
  fieldsError::Union{Array{typeof(1.0u"T"), 2}, Nothing} = nothing
  "Frequency of the measured field."
  fieldsFrequency::Union{Vector{typeof(1.0u"Hz")}, Nothing} = nothing
  "Applied current while measuring the matching position."
  currents::Union{Array{typeof(1.0u"A"), 2}, Nothing} = nothing
  "Start time of the measurement or, if defined as a vector, the timestamp of the matching position."
  timestamp::Union{DateTime, Vector{DateTime}, Nothing} = nothing
  "Offset of the hall probe's active areas."
  sensorOffset::Union{Vector{typeof(1.0u"m")}, Nothing} = nothing
  "Temperature at the start of the measurement or, if defined as a vector, the temperature when measuring the matching position."
  temperature::Union{typeof(1.0u"°C"), Vector{typeof(1.0u"°C")}, Nothing} = nothing
end

function MagneticFieldMeasurement(filename::String)
  return loadMagneticFieldMeasurement(filename)
end

function MagneticFieldMeasurement(file::HDF5.File)
  return loadMagneticFieldMeasurement(file)
end

function saveMagneticFieldAsHDF5(measurement::MagneticFieldMeasurement, filename::String)
  h5open(filename, "w") do file
    saveMagneticFieldAsHDF5(measurement, file)
  end
end

function saveMagneticFieldAsHDF5(measurement::MagneticFieldMeasurement, file::HDF5.File)
  write(file, measurement.positions)

  if !isnothing(description(measurement))
    write(file, "/description", description(measurement))
  end

  if !ismissing(fields(measurement))
    write(file, "/fields", ustrip.(fields(measurement)))
  end

  if !isnothing(fieldsError(measurement))
    write(file, "/fieldsError", ustrip.(fieldsError(measurement)))
  end

  if !isnothing(fieldsFrequency(measurement))
    write(file, "/fieldsFrequency", ustrip.(fieldsFrequency(measurement)))
  end

  if !isnothing(currents(measurement))
    write(file, "/currents", ustrip.(currents(measurement)))
  end

  if !isnothing(timestamp(measurement))
    write(file, "/timestamp", string.(timestamp(measurement)))
  end

  if !isnothing(sensorOffset(measurement))
    write(file, "/sensorOffset", ustrip.(sensorOffset(measurement)))
  end

  if !isnothing(temperature(measurement))
    write(file, "/temperature", ustrip.(temperature(measurement)))
  end
end

function loadMagneticFieldMeasurement(filename::String)
  h5open(filename, "r") do file
    return loadMagneticFieldMeasurement(file)
  end
end

function loadMagneticFieldMeasurement(file::HDF5.File)
  splattingDict = Dict{Symbol, Any}()
  splattingDict[:positions] = Positions(file)

  if haskey(file, "fields")
    splattingDict[:fields] = read(file, "fields")u"T"
  else
    error("The HDF5 file for a magnetic field measurement must contain a field vector.")
  end

  if typeof(splattingDict[:positions]) == MeanderingGridPositions
    splattingDict[:fields] = splattingDict[:fields][:, getPermutation(splattingDict[:positions]), :]
    splattingDict[:positions] = splattingDict[:positions].grid
  end

  if haskey(file, "description")
    splattingDict[:description] = read(file, "description")
  end

  if haskey(file, "fieldsError")
    splattingDict[:fieldsError] = read(file, "fieldsError")u"T"
  end

  if haskey(file, "fieldsFrequency")
    splattingDict[:fieldsFrequency] = read(file, "fieldsFrequency")u"Hz"
  end

  if haskey(file, "currents")
    splattingDict[:currents] = read(file, "currents")u"A"
  end

  if haskey(file, "timestamp")
    splattingDict[:timestamp] = DateTime(read(file, "timestamp"))
  end

  if haskey(file, "sensorOffset")
    splattingDict[:sensorOffset] = read(file, "sensorOffset")u"m"
  end

  if haskey(file, "temperature")
    splattingDict[:temperature] = read(file, "temperature")u"°C"
  end

  return MagneticFieldMeasurement(;splattingDict...)
end

# Alias for backwards compatability
loadMagneticField(filename::String) = loadMagneticFieldMeasurement(filename)

# Create getter and setter for all fields of `MagneticFieldMeasurement`
for (fieldname, fieldtype) in zip(fieldnames(MagneticFieldMeasurement), fieldtypes(MagneticFieldMeasurement))
  fieldnameStr = string(fieldname)

  # At the moment, this should be a Union
  missingOrNothing = (fieldtype.b <: Union{Missing, Nothing}) ? fieldtype.b : fieldtype.a
  fieldtype = (fieldtype.b <: Union{Missing, Nothing}) ? fieldtype.a : fieldtype.b

  #@info "" fieldnameStr missingOrNothing fieldtype

  @eval begin
    export $fieldname

    # TODO: Add docstring from struct; I did not yet find a way to retrieve it
    function $(fieldname)(measurement::MagneticFieldMeasurement)::Union{$fieldtype, $missingOrNothing}
      return measurement.$fieldname
    end
  
    # TODO: Add docstring from struct; I did not yet find a way to retrieve it
    function $(fieldname)(measurement::MagneticFieldMeasurement, value::Union{$fieldtype, $missingOrNothing})
      measurement.$fieldname = value
    end
  end
end