module MPIFiles

using UUIDs
using Graphics: @mustimplement

using LinearOperatorCollection
using FFTW
using AxisArrays
const axes = Base.axes
using Interpolations
using HDF5
using Dates
using DelimitedFiles
using ImageMetadata
using ImageAxes
using LinearAlgebra
using Random
using Mmap
using Statistics
using Unitful
using CodecZlib
using Tar
using Pkg.PlatformEngines
using Pkg.GitTools
using Pkg.Artifacts
using Unitful
using InteractiveUtils
using UnitfulAngles
using UnitfulParsableString
using Inflate, SHA
using StableRNGs
using REPL: fielddoc
using DocStringExtensions

### global import list ###

import Base: convert, get, getindex, haskey, iterate, length, ndims, range, read, show, time, write, close
import FileIO: save
import HDF5: h5read
import Interpolations: interpolate

### export list ###

export MPIFile

# general parameters
export version, uuid

# study parameters
export studyName, studyNumber, studyUuid, studyDescription, studyTime

# experiment parameters
export experimentName, experimentNumber, experimentUuid, experimentDescription, experimentSubject,
  experimentIsSimulation, experimentIsCalibration,
  experimentHasMeasurement, experimentHasReconstruction

# tracer parameters
export tracerName, tracerBatch, tracerVolume, tracerConcentration,
  tracerSolute, tracerInjectionTime, tracerVendor

# scanner parameters
export scannerBoreSize, scannerFacility, scannerOperator, scannerManufacturer,
   scannerName, scannerTopology

# acquisition parameters
export acqStartTime, acqNumFrames, acqNumAverages,
  acqGradient, acqOffsetField, acqNumPeriodsPerFrame

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
  dfDivider, dfWaveform, dfCycle

# receiver parameters
export rxNumChannels, rxBandwidth, rxNumSamplingPoints,
  rxTransferFunction, rxTransferFunctionFileName, rxHasTransferFunction, rxUnit,
  rxDataConversionFactor, rxInductionFactor

# measurements
export measData, measDataTDPeriods, measIsFourierTransformed, measIsTFCorrected,
  measIsTransferFunctionCorrected,
  measIsBGCorrected, measIsBackgroundCorrected, measIsFastFrameAxis,
  measIsFramePermutation, measIsFrequencySelection,
  measIsBGFrame, measIsBackgroundFrame, measIsSpectralLeakageCorrected, measFramePermutation,
  measFrequencySelection, measIsSparsityTransformed, measIsCalibProcessed

# calibrations
export calibSNR, calibSnr, calibFov, calibFieldOfView, calibFovCenter,
  calibFieldOfViewCenter, calibSize, calibOrder, calibPositions,
  calibOffsetFields, calibDeltaSampleSize,
  calibMethod, calibIsMeanderingGrid

# reconstruction results
export recoData, recoFov, recoFieldOfView, recoFovCenter, recoFieldOfViewCenter,
  recoSize, recoOrder, recoPositions

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
@mustimplement studyTime(f::MPIFile)

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
@mustimplement scannerBoreSize(f::MPIFile)
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
@mustimplement rxTransferFunctionFileName(f::MPIFile)
@mustimplement rxHasTransferFunction(f::MPIFile)
@mustimplement rxInductionFactor(f::MPIFile)
@mustimplement rxUnit(f::MPIFile)
@mustimplement rxDataConversionFactor(f::MPIFile)

# measurements
@mustimplement measData(f::MPIFile)
@mustimplement measDataTD(f::MPIFile)
@mustimplement measDataTDPeriods(f::MPIFile, periods)
@mustimplement measFramePermutation(f::MPIFile)
@mustimplement measFrequencySelection(f::MPIFile)
@mustimplement measIsBGCorrected(f::MPIFile)
@mustimplement measIsBGFrame(f::MPIFile)
@mustimplement measIsFastFrameAxis(f::MPIFile)
@mustimplement measIsFourierTransformed(f::MPIFile)
@mustimplement measIsFramePermutation(f::MPIFile)
@mustimplement measIsFrequencySelection(f::MPIFile)
@mustimplement measIsSparsityTransformed(f::MPIFile)
@mustimplement measIsSpectralLeakageCorrected(f::MPIFile)
@mustimplement measIsTFCorrected(f::MPIFile)
@mustimplement measSparsityTransformation(f::MPIFile)
@mustimplement measSubsamplingIndices(f::MPIFile)
@mustimplement measIsCalibProcessed(b::MPIFile)

# calibrations
@mustimplement calibSNR(f::MPIFile)
@mustimplement calibFov(f::MPIFile)
@mustimplement calibFovCenter(f::MPIFile)
@mustimplement calibSize(f::MPIFile)
@mustimplement calibOrder(f::MPIFile)
@mustimplement calibPositions(f::MPIFile)
@mustimplement calibOffsetFields(f::MPIFile)
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
@mustimplement recoIsOverscanRegion(f::MPIFile)

# additional functions that should be implemented by an MPIFile
@mustimplement filepath(f::MPIFile)


include("Derived.jl")
include("Custom.jl")
include("Utils.jl")
include("FramePermutation.jl")

### Concrete implementations ###
include("MDF.jl")
include("MDFInMemory.jl")
include("MDFCommon.jl")
include("Brukerfile.jl")
include("IMT.jl")

# This dispatches on the file extension and automatically
# generates the correct type
function MPIFile(filename::AbstractString; kargs...)
  filenamebase, ext = splitext(filename)
  if ext == ".mdf" || ext == ".hdf" || ext == ".h5"
    file = h5open(filename, "r")
    if haskey(file, "/version")
      return MDFFile(filename, file) # MDFFile currently has no kargs
    else
      return IMTFile(filename, file; kargs...)
    end
  else
    if isfile(joinpath(filename, "mdf"))
      filenameMDF = readline(joinpath(filename, "mdf"))
      return MDFFile(filenameMDF)
    else
      return BrukerFile(filename; kargs...)
    end
  end
end

export DMPIFile
"""
    DMPIFile(args...; worker)
  
Construct a distributed `MPIFile` using the given `args` on the given `worker` process.

See `MPIFile`
"""
function DMPIFile end

function show(io::IO, f::MPIFile)
  print(io, supertype(typeof(f)))
  print(io, "\n\tStudy: ")
  show(io, studyName(f))
  print(io, ", ")
  show(io, studyTime(f))
  print(io, "\n\tExperiment: ")
  show(io, experimentName(f))
  print(io, ", ")
  show(io, acqStartTime(f))
  print(io, "\n")
end

# Opens a set of MPIFiles
function MPIFile(filenames::Vector)
  return map(x -> MPIFile(x), filenames)
end

# For the do block
function MPIFile(h::Function, args...; kargs...)
  f = MPIFile(args...; kargs...)
  try
    h(f)
  finally
    close(f)
  end
end

Base.length(f::MPIFile) = 1
Base.close(f::MPIFile) = nothing

include("TransferFunction.jl")
include("MultiMPIFile.jl")
include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Conversion.jl")
include("RecoData.jl")
include("DatasetStore/DatasetStore.jl")
include("MixingFactors.jl")
include("Positions/Positions.jl")
include("MagneticFieldMeasurement.jl")

end # module
