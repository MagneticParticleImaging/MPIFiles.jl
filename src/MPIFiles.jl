__precompile__()
module MPIFiles

using Compat
using ProgressMeter
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
      experimentIsSimulation, experimentIsCalibration,
      experimentHasMeasurement, experimentHasReconstruction

# tracer parameters
export tracerName, tracerBatch, tracerVolume, tracerConcentration,
       tracerSolute, tracerInjectionTime, tracerVendor

# scanner parameters
export scannerFacility, scannerOperator, scannerManufacturer, scannerName,
       scannerTopology

# acquisition parameters
export acqStartTime, acqFramePeriod, acqNumFrames, acqNumAverages,
       acqGradient, acqOffsetField, acqOffsetFieldShift, acqNumPeriods, acqSize

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfPeriod

# receiver parameters
export rxNumChannels, rxBandwidth, rxNumSamplingPoints,
       rxTransferFunction, rxUnit, rxDataConversionFactor, rxInductionFactor

# measurements
export measData, measIsFourierTransformed, measIsTFCorrected,
       measIsBGCorrected, measIsTransposed,
       measIsFramePermutation, measIsFrequencySelection,
       measIsBGFrame, measIsSpectralLeakageCorrected, measFramePermutation

# calibrations
export calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

# additional functions that should be implemented by an MPIFile
export filepath, systemMatrixWithBG, systemMatrix

export selectedChannels
### Interface of an MPIFile ###

@compat abstract type MPIFile end

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
@mustimplement scannerName(f::MPIFile)
@mustimplement scannerTopology(f::MPIFile)

# acquisition parameters
@mustimplement acqStartTime(f::MPIFile)
@mustimplement acqFramePeriod(f::MPIFile)
@mustimplement acqNumAverages(f::MPIFile)
@mustimplement acqNumPeriods(f::MPIFile)
@mustimplement acqNumFrames(f::MPIFile)
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
@mustimplement rxBandwidth(f::MPIFile)
@mustimplement rxNumSamplingPoints(f::MPIFile)
@mustimplement rxTransferFunction(f::MPIFile)
@mustimplement rxInductionFactor(f::MPIFile)
@mustimplement rxUnit(f::MPIFile)
@mustimplement rxDataConversionFactor(f::MPIFile)

# measurements
@mustimplement measData(f::MPIFile)
@mustimplement measIsSpectralLeakageCorrected(f::MPIFile)
@mustimplement measIsFourierTransformed(f::MPIFile)
@mustimplement measIsTFCorrected(f::MPIFile)
@mustimplement measIsFrequencySelecton(f::MPIFile)
@mustimplement measIsBGCorrected(f::MPIFile)
@mustimplement measIsTransposed(f::MPIFile)
@mustimplement measIsFramePermutation(f::MPIFile)
@mustimplement measIsBGFrame(f::MPIFile)
@mustimplement measFramePermutation(f::MPIFile)

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
export rxNumFrequencies, acqFov, acqFovCenter, rxFrequencies, rxTimePoints
rxNumFrequencies(f::MPIFile) = floor(Int,rxNumSamplingPoints(f) ./ 2 .+ 1)
function rxFrequencies(f::MPIFile)
  numFreq = rxNumFrequencies(f)
  a = collect(0:(numFreq-1))./(numFreq-1).*rxBandwidth(f)
  return a
end
function rxTimePoints(f::MPIFile)
  numTP = rxNumSamplingPoints(f)
  a = collect(0:(numTP-1))./(numTP).*dfPeriod(f)
  return a
end
function acqFov(f::MPIFile)
 return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradient(f) )
end
function acqFovCenter(f::MPIFile)
 return acqOffsetField(f) ./ abs.( acqGradient(f) )
end

export acqNumFGFrames, acqNumBGFrames, measFGFrameIdx, measBGFrameIdx

acqNumFGFrames(f::MPIFile) = acqNumFrames(f) - acqNumBGFrames(f)
acqNumBGFrames(f::MPIFile) = sum(measIsBGFrame(f))

function measBGFrameIdx(f::MPIFile)
  idx = zeros(Int64, acqNumBGFrames(f))
  j = 1
  mask = measIsBGFrame(f)
  for i=1:acqNumFrames(f)
    if mask[i]
      idx[j] = i
      j += 1
    end
  end
  return idx
end

function measFGFrameIdx(f::MPIFile)
  mask = measIsBGFrame(f)
  if !any(mask)
    #shortcut
    return 1:acqNumFrames(f)
  end
  idx = zeros(Int64, acqNumFGFrames(f))
  j = 1
  for i=1:acqNumFrames(f)
    if !mask[i]
      idx[j] = i
      j += 1
    end
  end
  return idx
end

# We assume that systemMatrixWithBG has already reordered the BG data
# to the end
systemMatrix(f::MPIFile) = systemMatrixWithBG(f)[1:acqNumFGFrames(f),:,:,:]

### Concrete implementations ###

include("Custom.jl")
include("MDF.jl")
include("Brukerfile.jl")

# This dispatches on the file extension and automatically
# generates the correct type
function (::Type{MPIFile})(filename::AbstractString; kargs...)
  filenamebase, ext = splitext(filename)
  if ext == ".mdf" || ext == ".hdf" || ext == ".h5"
    return MDFFile(filename; kargs...)
  else
    return BrukerFile(filename; kargs...)
  end
end

# Opens a set of MPIFiles
function (::Type{MPIFile})(filenames::Vector)
  return map(x->MPIFile(x),filenames)
end

optParam(param, default) = (param == nothing) ? default : param

# Support for handling complex datatypes in HDF5 files
function writeComplexArray{T,D}(file, dataset, A::Array{Complex{T},D})
  d_type_compound = HDF5.h5t_create(HDF5.H5T_COMPOUND,2*sizeof(T))
  HDF5.h5t_insert(d_type_compound, "r", 0 , HDF5.hdf5_type_id(T))
  HDF5.h5t_insert(d_type_compound, "i", sizeof(T) , HDF5.hdf5_type_id(T))

  shape = collect(reverse(size(A)))
  space = HDF5.h5s_create_simple(D, shape, shape)

  dset_compound = HDF5.h5d_create(file, dataset, d_type_compound, space,
                                  HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT)
  HDF5.h5s_close(space)

  HDF5.h5d_write(dset_compound, d_type_compound, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, A)

  HDF5.h5d_close(dset_compound)
  HDF5.h5t_close(d_type_compound)
end

function isComplexArray(file, dataset)
  if eltype(file[dataset]) <: HDF5.HDF5Compound{2}
    if HDF5.h5t_get_member_name(datatype(file[dataset]).id,0) == "r" &&
      HDF5.h5t_get_member_name(datatype(file[dataset]).id,1) == "i"
        return true
    end
  end
  return false
end

function getComplexType(file, dataset)
  T = HDF5.hdf5_to_julia_eltype(
            HDF5Datatype(
              HDF5.h5t_get_member_type( datatype(file[dataset]).id, 0 )
          )
        )
    return Complex{T}
end

function readComplexArray(file::HDF5File, dataset)
  T = getComplexType(file, dataset)
  A = copy(readmmap(file[dataset],Array{getComplexType(file,dataset)}))
  return A
end

function readComplexArray(filename::String, dataset)
  h5open(filename, "r") do file
    return readComplexArray(file, dataset)
  end
end

include("MultiMPIFile.jl")
include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Conversion.jl")
include("Image.jl")
include("DatasetStore.jl")
include("Positions.jl")

end # module
