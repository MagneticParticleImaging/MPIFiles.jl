module MPIFiles

using Graphics: @mustimplement

import Base: ndims, time, show, getindex

### export list ###

export MPIFile

# general parameters
export version, uuid

# study parameters
export studyName, studyNumber, studyUuid, studyDescription

# experiment parameters
export experimentName, experimentNumber, experimentUuid, experimentDescription, experimentSubject,
      experimentIsSimulation, experimentIsCalibration, experimentHasProcessing

# tracer parameters
export tracerName, tracerBatch, tracerVolume, tracerConcentration,
       tracerSolute, tracerInjectionTime, tracerVendor

# scanner parameters
export scannerFacility, scannerOperator, scannerManufacturer, scannerModel,
       scannerTopology

# acquisition parameters
export acqStartTime, acqNumFrames, acqNumBGFrames, acqFramePeriod, acqNumPatches,
       acqGradient, acqOffsetField, acqOffsetFieldShift

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfPeriod

# receiver parameters
export rxNumChannels, rxNumAverages, rxBandwidth, rxNumSamplingPoints,
       rxTransferFunction

# measurements
export measUnit, measDataConversionFactor, measData, measBGData, 
       measDataTimeOrder, measBGDataTimeOrder

# processing
export procData, procIsFourierTransformed, procIsTFCorrected, procIsAveraged,
       procIsFramesSelected, procIsBGCorrected, procIsTransposed,
       procFramePermutation

# calibrations
export calibSystemMatrixData, calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

# additional functions that should be implemented by an MPIFile
export filepath


### Interface of an MPIFile ###

abstract MPIFile

# general parameters
@mustimplement version(f::MPIFile)
@mustimplement uuid(f::MPIFile)
@mustimplement time(f::MPIFile)

# study parameters
@mustimplement studyName(f::MPIFile)
@mustimplement studyNumber(f::MPIFile)
@mustimplement studyUuid(f::MPIFile)
@mustimplement studyDescription(f::MPIFile)

# experiment parameters
@mustimplement experimentName(f::MPIFile)
@mustimplement experimentNumber(f::MPIFile)
@mustimplement experimentUuid(f::MPIFile)
@mustimplement experimentDescription(f::MPIFile)
@mustimplement experimentSubject(f::MPIFile)
@mustimplement experimentIsSimulation(f::MPIFile)
@mustimplement experimentIsCalibration(f::MPIFile)
@mustimplement experimentHasProcessing(f::MPIFile)
@mustimplement experimentHasReconstruction(f::MPIFile)
@mustimplement experimentHasMeasurement(f::MPIFile)

# tracer parameters
@mustimplement tracerName(f::MPIFile)
@mustimplement tracerBatch(f::MPIFile)
@mustimplement tracerVolume(f::MPIFile)
@mustimplement tracerConcentration(f::MPIFile)
@mustimplement tracerSolute(f::MPIFile)
@mustimplement tracerInjectionTime(f::MPIFile)

# scanner parameters
@mustimplement scannerFacility(f::MPIFile)
@mustimplement scannerOperator(f::MPIFile)
@mustimplement scannerManufacturer(f::MPIFile)
@mustimplement scannerModel(f::MPIFile)
@mustimplement scannerTopology(f::MPIFile)

# acquisition parameters
@mustimplement acqStartTime(f::MPIFile)
@mustimplement acqNumFrames(f::MPIFile)
@mustimplement acqNumBGFrames(f::MPIFile)
@mustimplement acqFramePeriod(f::MPIFile)
@mustimplement acqNumPatches(f::MPIFile)
@mustimplement acqGradient(f::MPIFile)
@mustimplement acqOffsetField(f::MPIFile)
@mustimplement acqOffsetFieldShift(f::MPIFile)

# drive-field parameters
@mustimplement dfNumChannels(f::MPIFile)
@mustimplement dfStrength(f::MPIFile)
@mustimplement dfPhase(f::MPIFile)
@mustimplement dfBaseFrequency(f::MPIFile)
@mustimplement dfCustomWaveform(f::MPIFile)
@mustimplement dfDivider(f::MPIFile)
@mustimplement dfWaveform(f::MPIFile)
@mustimplement dfPeriod(f::MPIFile)

# receiver properties
@mustimplement rxNumChannels(f::MPIFile)
@mustimplement rxNumAverages(f::MPIFile)
@mustimplement rxBandwidth(f::MPIFile)
@mustimplement rxNumSamplingPoints(f::MPIFile)
@mustimplement rxTransferFunction(f::MPIFile)

# measurements
@mustimplement measUnit(f::MPIFile)
@mustimplement measDataConversionFactor(f::MPIFile)
@mustimplement measData(f::MPIFile)
@mustimplement measBGData(f::MPIFile)
@mustimplement measDataTimeOrder(f::MPIFile)
@mustimplement measBGDataTimeOrder(f::MPIFile)

# processing
@mustimplement procData(f::MPIFile)
@mustimplement procData(f::MPIFile, frequencies)
@mustimplement procIsFourierTransformed(f::MPIFile)
@mustimplement procIsTFCorrected(f::MPIFile)
@mustimplement procIsAveraged(f::MPIFile)
@mustimplement procIsFramesSelected(f::MPIFile)
@mustimplement procIsBGCorrected(f::MPIFile)
@mustimplement procIsTransposed(f::MPIFile)
@mustimplement procFramePermutation(f::MPIFile)

# calibrations
@mustimplement calibSNR(f::MPIFile)
@mustimplement calibFov(f::MPIFile)
@mustimplement calibFovCenter(f::MPIFile)
@mustimplement calibSize(f::MPIFile)
@mustimplement calibOrder(f::MPIFile)
@mustimplement calibPositions(f::MPIFile)
@mustimplement calibOffsetField(f::MPIFile)
@mustimplement calibDeltaSampleSize(f::MPIFile)
@mustimplement calibMethod(f::MPIFile)

# reconstruction results
@mustimplement recoData(f::MPIFile)
@mustimplement recoFov(f::MPIFile)
@mustimplement recoFovCenter(f::MPIFile)
@mustimplement recoSize(f::MPIFile)
@mustimplement recoOrder(f::MPIFile)
@mustimplement recoPositions(f::MPIFile)

# additional functions that should be implemented by an MPIFile
@mustimplement filepath(f::MPIFile)



function str2uuid(str::String)
  if contains(str,"-")
    str_ = str
  else
    str_ = string(str[1:8],"-",str[9:12],"-",str[13:16],"-",str[17:20],"-",str[21:end])
  end
  return Base.Random.UUID(str_)
end
str2uuid(str::Void) = str

#TODO Move to misc
export rxNumFrequencies, acqFov, acqFovCenter, rxFrequencies
rxNumFrequencies(f::MPIFile) = floor(Int,rxNumSamplingPoints(f) ./ 2 .+ 1)
function rxFrequencies(f::MPIFile)
  numFreq = rxNumFrequencies(f)
  a = collect(0:(numFreq-1))./(numFreq-1).*rxBandwidth(f)
  return a
end
function acqFov(f::MPIFile)
 return addLeadingSingleton( 2*vec(dfStrength(f)) ./ vec(abs( acqGradient(f) )),2)
end
function acqFovCenter(f::MPIFile)
 return addLeadingSingleton( vec(acqOffsetField(f)) ./ vec(abs( acqGradient(f) )),2)
end

function measData(f::MPIFile, args...; backgroundData=false, kargs...) 
  if backgroundData
    measBGData(f,args...;kargs...)
  else
    measData(f,args...;kargs...)
  end
end

function measDataConv(f::MPIFile, args...;backgroundData=false)
  data = measData(f, args..., backgroundData=backgroundData)
  a = measDataConversionFactor(f)
  data = map(Float32, data)
  scale!(data, a[1])
  data[:] .+= a[2]
  return data
end

### Concrete implementations ###

include("Custom.jl")
include("MDF.jl")
include("RawFile.jl")
include("Brukerfile.jl")

# This dispatches on the file extension and automatically
# generates the correct type
function (::Type{MPIFile})(filename::AbstractString)
  filenamebase, ext = splitext(filename)
  if ext == ".mdf" || ext == ".hdf" || ext == ".h5"
    return MDFFile(filename)
  else
    return BrukerFile(filename)
  end
end

# Opens a set of MPIFiles
function (::Type{MPIFile})(filenames::Vector)
  return map(x->MPIFile(x),filenames)
end

optParam(param, default) = (param == nothing) ? default : param

include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Conversion.jl")
include("Image.jl")

### Misc functions ###
#include("Misc.jl")


end # module
