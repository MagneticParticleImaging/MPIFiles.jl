module MPIFiles

using Graphics: @mustimplement

import Base: ndims, time, show, getindex

### export list ###

export MPIFile

# general parameters
export version, uuid

# study parameters
export studyName, studyExperiment, studyDescription, studySubject,
       studyIsSimulation, studyIsCalibration

# tracer parameters
export tracerName, tracerBatch, tracerVolume, tracerConcentration,
       tracerSolute, tracerInjectionTime, tracerVendor

# scanner parameters
export scannerFacility, scannerOperator, scannerManufacturer, scannerModel,
       scannerTopology

# acquisition parameters
export acqStartTime, acqNumFrames, acqNumBGFrames, acqFramePeriod, acqNumPatches,
       acqGradient, acqOffsetField, acqFov, acqFovCenter

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfPeriod

# receiver parameters
export rxNumChannels, rxNumAverages, rxBandwidth, rxNumSamplingPoints, 
       rxTransferFunction

# measurements
export measUnit, measRawDataConversion,
       measData, measDataTimeOrder, measBGData, measBGDataTimeOrder

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
@mustimplement studyExperiment(f::MPIFile)
@mustimplement studyDescription(f::MPIFile)
@mustimplement studySubject(f::MPIFile)
@mustimplement studyIsSimulation(f::MPIFile)
@mustimplement studyIsCalibration(f::MPIFile)

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
@mustimplement acqFov(f::MPIFile)
@mustimplement acqFovCenter(f::MPIFile)

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
@mustimplement measRawDataConversion(f::MPIFile)
@mustimplement measData(f::MPIFile)
@mustimplement measDataTimeOrder(f::MPIFile)
@mustimplement measBGData(f::MPIFile)
@mustimplement measBGDataTimeOrder(f::MPIFile)

# calibrations
@mustimplement calibSystemMatrixData(f::MPIFile)
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





### Concrete implementations ###
include("MDF.jl")

#TODO Move to misc
rxNumFrequencies(f::MPIFile) = floor(Int,rxNumSamplingPoints(f) ./ 2 .+ 1)
function rxFrequencies(f::MPIFile)
  numFreq = rxNumFrequencies(f)
  a = collect(0:(numFreq-1))./(numFreq-1).*rxBandwidth(b)
  return a
end


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

include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Conversion.jl")

### Misc functions ###
#include("Misc.jl")


end # module
