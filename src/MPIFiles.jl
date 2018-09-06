module MPIFiles

using Reexport

using Compat
@reexport using LinearAlgebra
@reexport using Statistics
@reexport using FFTW
@reexport using Random
using Compat.UUIDs
using ProgressMeter
using Graphics: @mustimplement
@reexport using ImageAxes
using AxisArrays
const axes = Base.axes
@reexport using ImageMetadata
@reexport using Unitful
using Interpolations
@reexport using Mmap
@reexport using Dates
@reexport using DelimitedFiles


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
export acqStartTime, acqNumFrames, acqNumAverages,
       acqGradient, acqOffsetField, acqNumPeriodsPerFrame, acqSize

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfCycle

# receiver parameters
export rxNumChannels, rxBandwidth, rxNumSamplingPoints,
       rxTransferFunction, rxUnit, rxDataConversionFactor, rxInductionFactor

# measurements
export measData, measDataTDPeriods, measIsFourierTransformed, measIsTFCorrected,
       measIsBGCorrected, measIsTransposed,
       measIsFramePermutation, measIsFrequencySelection,
       measIsBGFrame, measIsSpectralLeakageCorrected, measFramePermutation

# calibrations
export calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod, calibIsMeanderingGrid

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

# additional functions that should be implemented by an MPIFile
export filepath, systemMatrixWithBG, systemMatrix

export selectedChannels
### Interface of an MPIFile ###

abstract type MPIFile end

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
@mustimplement acqNumAverages(f::MPIFile)
@mustimplement acqNumPeriodsPerFrame(f::MPIFile)
@mustimplement acqNumFrames(f::MPIFile)
@mustimplement acqGradient(f::MPIFile)
@mustimplement acqOffsetField(f::MPIFile)

# drive-field parameters
@mustimplement dfNumChannels(f::MPIFile)
@mustimplement dfStrength(f::MPIFile)
@mustimplement dfPhase(f::MPIFile)
@mustimplement dfBaseFrequency(f::MPIFile)
@mustimplement dfCustomWaveform(f::MPIFile)
@mustimplement dfDivider(f::MPIFile)
@mustimplement dfWaveform(f::MPIFile)
@mustimplement dfCycle(f::MPIFile)

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
@mustimplement measDataTD(f::MPIFile)
@mustimplement measDataTDPeriods(f::MPIFile, periods)
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
@mustimplement calibIsMeanderingGrid(f::MPIFile)

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
  if occursin("-", str)
    str_ = str
  else
    str_ = string(str[1:8],"-",str[9:12],"-",str[13:16],"-",str[17:20],"-",str[21:end])
  end
  try
    u = UUID(str_)
    return u
  catch
    println("could not convert to UUID. str_= $(str_)  str=$(str) ")
    u = uuid4()
    return u
  end
end
str2uuid(str::Nothing) = str

# TODO Move the following functions to misc

export rxNumFrequencies, acqFov, rxFrequencies, rxTimePoints, acqGradientDiag

rxNumFrequencies(f::MPIFile) = floor(Int,rxNumSamplingPoints(f) ./ 2 .+ 1)

function rxFrequencies(f::MPIFile)
  numFreq = rxNumFrequencies(f)
  a = collect(0:(numFreq-1))./(numFreq-1).*rxBandwidth(f)
  return a
end

function rxTimePoints(f::MPIFile)
  numTP = rxNumSamplingPoints(f)
  a = collect(0:(numTP-1))./(numTP).*dfCycle(f)
  return a
end

function acqGradientDiag(f::MPIFile)
  g = acqGradient(f)
  g_ = reshape(g,9,size(g,3),size(g,4))
  return g_[[1,5,9],:,:]
end

function acqFov(f::MPIFile)
  if size(dfStrength(f)[1,:,:],1) == 3
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[:,1,:] )
  else
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[1,1,1] )
  end
end

#function acqFovCenter(f::MPIFile)
# return acqOffsetField(f) ./ abs.( acqGradient(f) ) # why was the absolute value taken here?
#end

export acqNumFGFrames, acqNumBGFrames, measFGFrameIdx, measBGFrameIdx, acqOffsetFieldShift,
       acqFramePeriod, acqNumPeriods, acqNumPatches, acqNumPeriodsPerPatch, measBGFrameBlockLengths

acqFramePeriod(b::MPIFile) = dfCycle(b) * acqNumAverages(b) * acqNumPeriodsPerFrame(b)

# numPeriods is the total number of DF periods in a measurement.
acqNumPeriods(f::MPIFile) = acqNumFrames(f)*acqNumPeriodsPerFrame(f)

function acqOffsetFieldShift(f::MPIFile)
    return acqOffsetField(f) ./ reshape( acqGradient(f),9,1,:)[[1,5,9],:,:]
end

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


function measBGFrameBlockLengths(mask)
  len = Vector{Int}(undef,0)

  groupIdxStart = -1
  for i=1:(length(mask)+1)
    if i <= length(mask) && mask[i] && groupIdxStart == -1
      groupIdxStart = i
    end
    if groupIdxStart != -1 && ((i == length(mask)+1) || !mask[i])
      push!(len, i-groupIdxStart)
      groupIdxStart = - 1
    end
  end
  return len
end


function acqNumPatches(f::MPIFile)
  # not valid for varying gradients / multi gradient
  shifts = acqOffsetFieldShift(f)
  return size(unique(shifts,dims=2),2)
end

function acqNumPeriodsPerPatch(f::MPIFile)
  return div(acqNumPeriodsPerFrame(f), acqNumPatches(f))
end

export unflattenOffsetFieldShift

unflattenOffsetFieldShift(f::MPIFile) = analyseFFPos(acqOffsetFieldShift(f))
function unflattenOffsetFieldShift(shifts::Array)
  # not valid for varying gradients / multi gradient
  uniqueShifts = unique(shifts, dims=2)
  numPeriodsPerFrame = size(shifts,2)
  numPatches = size(uniqueShifts,2)
  numPeriodsPerPatch = div(numPeriodsPerFrame, numPatches)

  allPeriods = 1:numPeriodsPerFrame

  flatIndices = zeros(Int64,numPatches,numPeriodsPerPatch)

  for i=1:numPatches
    flatIndices[i,:] = allPeriods[vec(sum(shifts .== uniqueShifts[:,i],dims=1)).==3]
  end

  return flatIndices
end

# We assume that systemMatrixWithBG has already reordered the BG data
# to the end
systemMatrix(f::MPIFile) = systemMatrixWithBG(f)[1:acqNumFGFrames(f),:,:,:]

function measDataTD(f, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))

  data1 = measData(f,frames,periods,receivers)

  if measIsTransposed(f)
    data2 = permutedims(data1, invperm([4,1,2,3]))
  else
    data2 = data1
  end

  if measIsFourierTransformed(f)
    dataTD = irfft(data2,2*size(data2,1)-1, 1)
  else
    dataTD = data2
  end
  return dataTD
end


### Concrete implementations ###

include("Custom.jl")
include("FramePermutation.jl")
include("MDF.jl")
include("Brukerfile.jl")

# This dispatches on the file extension and automatically
# generates the correct type
function MPIFile(filename::AbstractString; kargs...)
  filenamebase, ext = splitext(filename)
  if ext == ".mdf" || ext == ".hdf" || ext == ".h5"
    return MDFFile(filename; kargs...)
  else
    return BrukerFile(filename; kargs...)
  end
end

# Opens a set of MPIFiles
function MPIFile(filenames::Vector)
  return map(x->MPIFile(x),filenames)
end

# Support for handling complex datatypes in HDF5 files
function writeComplexArray(file, dataset, A::AbstractArray{Complex{T},D}) where {T,D}
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
include("MultiPatchMPIFile.jl")
include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Conversion.jl")
include("Image.jl")
include("DatasetStore.jl")
include("MixingFactors.jl")
include("positions/Positions.jl")


end # module
