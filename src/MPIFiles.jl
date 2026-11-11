module MPIFiles

using UUIDs
using Graphics: @mustimplement
using RelocatableFolders

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
using StaticArrays
using REPL: fielddoc
using DocStringExtensions
using NearestNeighbors

### global import list ###

import Base: convert, get, getindex, haskey, iterate, length, ndims, range, read, show, time, write, close
import Base: eltype, parent, parentindices, view
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
  measFrequencySelection, measIsSparsityTransformed, measSparsityTransformation, measSubsamplingIndices,
  measIsCalibProcessed

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
"""
    version(f::MPIFile) -> VersionNumber

Version of the file format. MDF: `/version`.
"""
@mustimplement version(f::MPIFile)
"""
    uuid(f::MPIFile) -> UUID

Universally Unique Identifier (RFC 4122) of the file. MDF: `/uuid`.
"""
@mustimplement uuid(f::MPIFile)
"""
    time(f::MPIFile) -> DateTime

UTC creation time of the data set. MDF: `/time`.
"""
@mustimplement time(f::MPIFile)

# study parameters
"""
    studyName(f::MPIFile) -> String

Name of the study. MDF: `/study/name`.
"""
@mustimplement studyName(f::MPIFile)
"""
    studyNumber(f::MPIFile) -> Integer

Number of the study. MDF: `/study/number`.
"""
@mustimplement studyNumber(f::MPIFile)
"""
    studyUuid(f::MPIFile) -> UUID

Universally Unique Identifier (RFC 4122) of the study. MDF: `/study/uuid`.
"""
@mustimplement studyUuid(f::MPIFile)
"""
    studyDescription(f::MPIFile) -> String

Short description of the study. MDF: `/study/description`.
"""
@mustimplement studyDescription(f::MPIFile)
"""
    studyTime(f::MPIFile) -> DateTime

UTC creation time of the study. MDF: `/study/time`.
"""
@mustimplement studyTime(f::MPIFile)

# experiment parameters
"""
    experimentName(f::MPIFile) -> String

Name of the experiment. MDF: `/experiment/name`.
"""
@mustimplement experimentName(f::MPIFile)
"""
    experimentNumber(f::MPIFile) -> Integer

Number of the experiment within the study. MDF: `/experiment/number`.
"""
@mustimplement experimentNumber(f::MPIFile)
"""
    experimentUuid(f::MPIFile) -> UUID

Universally Unique Identifier (RFC 4122) of the experiment. MDF: `/experiment/uuid`.
"""
@mustimplement experimentUuid(f::MPIFile)
"""
    experimentDescription(f::MPIFile) -> String

Short description of the experiment. MDF: `/experiment/description`.
"""
@mustimplement experimentDescription(f::MPIFile)
"""
    experimentSubject(f::MPIFile) -> String

Name of the subject that was imaged. MDF: `/experiment/subject`.
"""
@mustimplement experimentSubject(f::MPIFile)
"""
    experimentIsSimulation(f::MPIFile) -> Bool

Flag indicating if the data in this file is simulated rather than measured. MDF: `/experiment/isSimulation`.
"""
@mustimplement experimentIsSimulation(f::MPIFile)
"""
    experimentIsCalibration(f::MPIFile) -> Bool

Flag indicating if the file contains a calibration measurement, i.e. a system matrix.
"""
@mustimplement experimentIsCalibration(f::MPIFile)
"""
    experimentHasReconstruction(f::MPIFile) -> Bool

Flag indicating if the file contains reconstruction results.
"""
@mustimplement experimentHasReconstruction(f::MPIFile)
"""
    experimentHasMeasurement(f::MPIFile) -> Bool

Flag indicating if the file contains measurement data.
"""
@mustimplement experimentHasMeasurement(f::MPIFile)

# tracer parameters
"""
    tracerName(f::MPIFile) -> Vector{String}

Name of the tracer used in the experiment, one entry per tracer (A). MDF: `/tracer/name`.
"""
@mustimplement tracerName(f::MPIFile)
"""
    tracerBatch(f::MPIFile) -> Vector{String}

Batch of the tracer, size (A). MDF: `/tracer/batch`.
"""
@mustimplement tracerBatch(f::MPIFile)
"""
    tracerVolume(f::MPIFile) -> Vector{Float64}

Total volume of applied tracer in L, size (A). MDF: `/tracer/volume`.
"""
@mustimplement tracerVolume(f::MPIFile)
"""
    tracerConcentration(f::MPIFile) -> Vector{Float64}

Molar concentration of the solute per litre in mol(solute)/L, size (A). MDF: `/tracer/concentration`.
"""
@mustimplement tracerConcentration(f::MPIFile)
"""
    tracerSolute(f::MPIFile) -> Vector{String}

Solute of the tracer, e.g. Fe, size (A). MDF: `/tracer/solute`.
"""
@mustimplement tracerSolute(f::MPIFile)
"""
    tracerInjectionTime(f::MPIFile) -> Vector{DateTime}

UTC time at which the tracer injection started, size (A). Optional. MDF: `/tracer/injectionTime`.
"""
@mustimplement tracerInjectionTime(f::MPIFile)

# scanner parameters
"""
    scannerBoreSize(f::MPIFile) -> Float64

Diameter of the bore in m. Optional. MDF: `/scanner/boreSize`.
"""
@mustimplement scannerBoreSize(f::MPIFile)
"""
    scannerFacility(f::MPIFile) -> String

Facility where the MPI scanner is installed. MDF: `/scanner/facility`.
"""
@mustimplement scannerFacility(f::MPIFile)
"""
    scannerOperator(f::MPIFile) -> String

User who operates the MPI scanner. MDF: `/scanner/operator`.
"""
@mustimplement scannerOperator(f::MPIFile)
"""
    scannerManufacturer(f::MPIFile) -> String

Scanner manufacturer. MDF: `/scanner/manufacturer`.
"""
@mustimplement scannerManufacturer(f::MPIFile)
"""
    scannerName(f::MPIFile) -> String

Scanner name. MDF: `/scanner/name`.
"""
@mustimplement scannerName(f::MPIFile)
"""
    scannerTopology(f::MPIFile) -> String

Scanner topology, e.g. FFP, FFL or MPS. MDF: `/scanner/topology`.
"""
@mustimplement scannerTopology(f::MPIFile)

# acquisition parameters
"""
    acqStartTime(f::MPIFile) -> DateTime

UTC start time of the MPI measurement. MDF: `/acquisition/startTime`.
"""
@mustimplement acqStartTime(f::MPIFile)
"""
    acqNumAverages(f::MPIFile) -> Integer

Number of block averages per drive-field period
"""
@mustimplement acqNumAverages(f::MPIFile)
"""
    acqNumPeriodsPerFrame(f::MPIFile) -> Integer

Number of drive-field periods within a frame, denoted by J. MDF: `/acquisition/numPeriodsPerFrame`.
"""
@mustimplement acqNumPeriodsPerFrame(f::MPIFile)
"""
    acqNumFrames(f::MPIFile) -> Integer

Number of available measurement frames, denoted by N. MDF: `/acquisition/numFrames`.
"""
@mustimplement acqNumFrames(f::MPIFile)
"""
    acqGradient(f::MPIFile) -> Array{Float64,4}

Gradient strength of the selection field in T m⁻¹ µ₀⁻¹, size (3, 3, Y, J). MDF: `/acquisition/gradient`.
"""
@mustimplement acqGradient(f::MPIFile)
"""
    acqOffsetField(f::MPIFile) -> Array{Float64,3}

Applied offset field in T µ₀⁻¹, size (3, Y, J). MDF: `/acquisition/offsetField`.
"""
@mustimplement acqOffsetField(f::MPIFile)

# drive-field parameters, group /acquisition/drivefield
"""
    dfNumChannels(f::MPIFile) -> Integer

Number of drive field channels, denoted by D
"""
@mustimplement dfNumChannels(f::MPIFile)
"""
    dfStrength(f::MPIFile) -> Array{Real}

Applied drive field strength in Tµ₀⁻¹, size (F, D, J) 
"""
@mustimplement dfStrength(f::MPIFile)
"""
    dfPhase(f::MPIFile) -> Array{Real}

Applied drive field phase ϕ in rad, size (F, D, J)
"""
@mustimplement dfPhase(f::MPIFile)
"""
    dfBaseFrequency(f::MPIFile) -> Real

Base frequency to derive drive field frequencies in Hz
"""
@mustimplement dfBaseFrequency(f::MPIFile)
"""
    dfCustomWaveform(f::MPIFile) -> String

Name of the custom waveform if [`dfWaveform`](@ref) is custom. Optional. MDF: `/acquisition/drivefield/customWaveform`.
"""
@mustimplement dfCustomWaveform(f::MPIFile)
"""
    dfDivider(f::MPIFile) -> Array{Integer}

Divider of the `dfBaseFrequency` to determine the drive field frequencies
"""
@mustimplement dfDivider(f::MPIFile)
"""
    dfWaveform(f::MPIFile) -> Array{String}

Waveform type: sine, triangle or custom, size (F, D)
"""
@mustimplement dfWaveform(f::MPIFile)
"""
    dfCycle(f::MPIFile) -> Real

Trajectory cycle is determined by `lcm(dfDivider)/dfBaseFrequency`. It will not change
when averaging was applied. The duration for measuring the V data points (i.e. the
drive-field period) is given by the product of `dfCycle` and `acqNumAverages`
"""
@mustimplement dfCycle(f::MPIFile)

# receiver properties
"""
    rxNumChannels(f::MPIFile) -> Integer

Number of receive channels, denoted by C. MDF: `/acquisition/receiver/numChannels`.
"""
@mustimplement rxNumChannels(f::MPIFile)
"""
    rxBandwidth(f::MPIFile) -> Float64

Bandwidth of the receiver unit in Hz. MDF: `/acquisition/receiver/bandwidth`.
"""
@mustimplement rxBandwidth(f::MPIFile)
"""
    rxNumSamplingPoints(f::MPIFile) -> Integer

Number of sampling points during one drive-field period, denoted by V. MDF: `/acquisition/receiver/numSamplingPoints`.
"""
@mustimplement rxNumSamplingPoints(f::MPIFile)
"""
    rxTransferFunction(f::MPIFile) -> Array{ComplexF64,2}

Transfer function of the receive channels in Fourier domain, size (K, C). `nothing` if the file contains no transfer function. MDF: `/acquisition/receiver/transferFunction`.
"""
@mustimplement rxTransferFunction(f::MPIFile)
"""
    rxTransferFunctionFileName(f::MPIFile) -> String

Name of the file the transfer function was loaded from. Non-standard, `nothing` if not available.
"""
@mustimplement rxTransferFunctionFileName(f::MPIFile)
"""
    rxHasTransferFunction(f::MPIFile) -> Bool

Flag indicating if the file contains a transfer function.
"""
@mustimplement rxHasTransferFunction(f::MPIFile)
"""
    rxInductionFactor(f::MPIFile) -> Vector{Float64}

Induction factor mapping the projection of the magnetic moment to the voltage in the receive coil, size (C). Optional. MDF: `/acquisition/receiver/inductionFactor`.
"""
@mustimplement rxInductionFactor(f::MPIFile)
"""
    rxUnit(f::MPIFile) -> String

SI unit of the measured quantity, usually V. MDF: `/acquisition/receiver/unit`.
"""
@mustimplement rxUnit(f::MPIFile)
"""
    rxDataConversionFactor(f::MPIFile) -> Array{Float64,2}

Dimensionless scaling factor and offset (a_c, b_c) to convert raw data into the physical quantity given by [`rxUnit`](@ref), size (2, C). Optional. MDF: `/acquisition/receiver/dataConversionFactor`.
"""
@mustimplement rxDataConversionFactor(f::MPIFile)

# measurements
"""
    measData(f::MPIFile) -> Array

Measured data at the processing stage described by the `measIs*` flags. For unprocessed data the size is (V, C, J, N). See [`getMeasurements`](@ref) and [`getMeasurementsFD`](@ref) for a high level interface. MDF: `/measurement/data`.
"""
@mustimplement measData(f::MPIFile)
"""
    measDataTD(f::MPIFile) -> Array

Measured data in time domain.
"""
@mustimplement measDataTD(f::MPIFile)
"""
    measDataTDPeriods(f::MPIFile, periods) -> Array

Measured time-domain data of the given `periods`, where the period and frame dimensions are merged.
"""
@mustimplement measDataTDPeriods(f::MPIFile, periods)
"""
    measFramePermutation(f::MPIFile) -> Vector{Int64}

Indices of the original frame order, size (N). Optional if `!measIsFramePermutation(f)`. MDF: `/measurement/framePermutation`.
"""
@mustimplement measFramePermutation(f::MPIFile)
"""
    measFrequencySelection(f::MPIFile) -> Vector{Int64}

Indices of the selected frequency components, size (K). Optional if `!measIsFrequencySelection(f)`. MDF: `/measurement/frequencySelection`.
"""
@mustimplement measFrequencySelection(f::MPIFile)
"""
    measIsBGCorrected(f::MPIFile) -> Bool

Flag, if the background has been subtracted. MDF: `/measurement/isBackgroundCorrected`.
"""
@mustimplement measIsBGCorrected(f::MPIFile)
"""
    measIsBGFrame(f::MPIFile) -> Vector{Bool}

Mask indicating for each of the N frames if it is a background measurement (`true`) or not. MDF: `/measurement/isBackgroundFrame`.
"""
@mustimplement measIsBGFrame(f::MPIFile)
"""
    measIsFastFrameAxis(f::MPIFile) -> Bool

Flag, if the frame dimension N has been moved to the first (fastest) dimension in Julia order. MDF: `/measurement/isFastFrameAxis`.
"""
@mustimplement measIsFastFrameAxis(f::MPIFile)
"""
    measIsFourierTransformed(f::MPIFile) -> Bool

Flag, if the data is stored in frequency space. MDF: `/measurement/isFourierTransformed`.
"""
@mustimplement measIsFourierTransformed(f::MPIFile)
"""
    measIsFramePermutation(f::MPIFile) -> Bool

Flag, if the order of frames has been changed, see [`measFramePermutation`](@ref). MDF: `/measurement/isFramePermutation`.
"""
@mustimplement measIsFramePermutation(f::MPIFile)
"""
    measIsFrequencySelection(f::MPIFile) -> Bool

Flag, if only a subset of frequencies has been selected and stored, see [`measFrequencySelection`](@ref). MDF: `/measurement/isFrequencySelection`.
"""
@mustimplement measIsFrequencySelection(f::MPIFile)
"""
    measIsSparsityTransformed(f::MPIFile) -> Bool

Flag, if the foreground frames are compressed along the frame dimension. MDF: `/measurement/isSparsityTransformed`.
"""
@mustimplement measIsSparsityTransformed(f::MPIFile)
"""
    measIsSpectralLeakageCorrected(f::MPIFile) -> Bool

Flag, if spectral leakage correction has been applied. MDF: `/measurement/isSpectralLeakageCorrected`.
"""
@mustimplement measIsSpectralLeakageCorrected(f::MPIFile)
"""
    measIsTFCorrected(f::MPIFile) -> Bool

Flag, if the data has been corrected by the transfer function, see [`rxTransferFunction`](@ref). MDF: `/measurement/isTransferFunctionCorrected`.
"""
@mustimplement measIsTFCorrected(f::MPIFile)
"""
    measSparsityTransformation(f::MPIFile) -> String

Name of the applied sparsity transformation. Optional if `!measIsSparsityTransformed(f)`. MDF: `/measurement/sparsityTransformation`.
"""
@mustimplement measSparsityTransformation(f::MPIFile)
"""
    measSubsamplingIndices(f::MPIFile) -> Array{Int64,4}

Subsampling indices of the sparsity transformation. Optional if `!measIsSparsityTransformed(f)`. MDF: `/measurement/subsamplingIndices`.
"""
@mustimplement measSubsamplingIndices(f::MPIFile)
"""
    measIsCalibProcessed(b::MPIFile) -> Bool

Flag, if the data is stored in the processed calibration layout, i.e. frame permuted, Fourier transformed and with fast frame axis.
"""
@mustimplement measIsCalibProcessed(b::MPIFile)

# calibrations
"""
    calibSNR(f::MPIFile) -> Array{Float64,3}

Signal-to-noise estimate of the recorded frequency components, size (K, C, J). Optional. MDF: `/calibration/snr`.
"""
@mustimplement calibSNR(f::MPIFile)
"""
    calibFov(f::MPIFile) -> Vector{Float64}

Field of view of the system matrix in m, size (3). Optional. MDF: `/calibration/fieldOfView`.
"""
@mustimplement calibFov(f::MPIFile)
"""
    calibFovCenter(f::MPIFile) -> Vector{Float64}

Center of the system matrix relative to the scanner origin in m, size (3). Optional. MDF: `/calibration/fieldOfViewCenter`.
"""
@mustimplement calibFovCenter(f::MPIFile)
"""
    calibSize(f::MPIFile) -> Vector{Int64}

Number of voxels in each dimension, the product is O. Optional. MDF: `/calibration/size`.
"""
@mustimplement calibSize(f::MPIFile)
"""
    calibOrder(f::MPIFile) -> String

Ordering of the dimensions, default is `"xyz"`. Optional. MDF: `/calibration/order`.
"""
@mustimplement calibOrder(f::MPIFile)
"""
    calibPositions(f::MPIFile) -> Array{Float64,2}

Position of each of the grid points in m, size (3, O). Optional. MDF: `/calibration/positions`.
"""
@mustimplement calibPositions(f::MPIFile)
"""
    calibOffsetFields(f::MPIFile) -> Array{Float64,2}

Applied offset field strength in T µ₀⁻¹ to emulate a spatial position, size (3, O). Optional. MDF: `/calibration/offsetFields`.
"""
@mustimplement calibOffsetFields(f::MPIFile)
"""
    calibDeltaSampleSize(f::MPIFile) -> Vector{Float64}

Size of the delta sample used for the calibration scan in m, size (3). Optional. MDF: `/calibration/deltaSampleSize`.
"""
@mustimplement calibDeltaSampleSize(f::MPIFile)
"""
    calibMethod(f::MPIFile) -> String

Method used to obtain the calibration data, e.g. `"robot"`, `"hybrid"` or `"simulation"`. MDF: `/calibration/method`.
"""
@mustimplement calibMethod(f::MPIFile)
"""
    calibIsMeanderingGrid(f::MPIFile) -> Bool

Flag, if the calibration grid is meandering. Optional. MDF: `/calibration/isMeanderingGrid`.
"""
@mustimplement calibIsMeanderingGrid(f::MPIFile)

# reconstruction results
"""
    recoData(f::MPIFile) -> Array

Reconstructed data, size (S, P, Q) with S channels, P voxels and Q frames. See [`loadRecoData`](@ref) for a high level interface. MDF: `/reconstruction/data`.
"""
@mustimplement recoData(f::MPIFile)
"""
    recoFov(f::MPIFile) -> Vector{Float64}

Field of view of the reconstructed data in m, size (3). Optional. MDF: `/reconstruction/fieldOfView`.
"""
@mustimplement recoFov(f::MPIFile)
"""
    recoFovCenter(f::MPIFile) -> Vector{Float64}

Center of the reconstructed data relative to the scanner origin in m, size (3). Optional. MDF: `/reconstruction/fieldOfViewCenter`.
"""
@mustimplement recoFovCenter(f::MPIFile)
"""
    recoSize(f::MPIFile) -> Vector{Int64}

Number of voxels in each dimension, the product is P. Optional. MDF: `/reconstruction/size`.
"""
@mustimplement recoSize(f::MPIFile)
"""
    recoOrder(f::MPIFile) -> String

Ordering of the dimensions, default is `"xyz"`. Optional. MDF: `/reconstruction/order`.
"""
@mustimplement recoOrder(f::MPIFile)
"""
    recoPositions(f::MPIFile) -> Array{Float64,2}

Position of each of the grid points in m, size (3, P). Optional. MDF: `/reconstruction/positions`.
"""
@mustimplement recoPositions(f::MPIFile)
"""
    recoIsOverscanRegion(f::MPIFile) -> Vector{Bool}

Mask indicating for each of the P voxels if it is part of the overscan region (`true`) or not. Optional. MDF: `/reconstruction/isOverscanRegion`.
"""
@mustimplement recoIsOverscanRegion(f::MPIFile)

# additional functions that should be implemented by an MPIFile
"""
    filepath(f::MPIFile) -> String

Path of the file or directory the data was loaded from.
"""
@mustimplement filepath(f::MPIFile)

abstract type AbstractValueDomain end
struct TimeDomain <: AbstractValueDomain end
struct FrequencyDomain <: AbstractValueDomain end

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
"""
    MPIFile(filename::AbstractString; kargs...) -> MPIFile
    MPIFile(filenames::Vector) -> Vector{<:MPIFile}
    MPIFile(f::Function, args...; kargs...)

Open an MPI data set and return a file handle of the matching type. Files with the extension
`.mdf`, `.hdf` or `.h5` are opened as MDF file (or IMT file if they contain no version), directories
are opened as Bruker data set. The data is loaded on demand when calling the accessor functions,
so opening a file is cheap.

The `do`-block form closes the file after `f` has been applied to it:
```julia
MPIFile("measurement.mdf") do f
  acqNumFrames(f)
end
```
"""
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

# Can't get https://github.com/JuliaLang/julia/pull/56368 to work, therefore we expor the abstract type here
# for other extensions/packages to dispatch on
export DMPIFile
"""
    DMPIFile(args...; worker)
  
Construct a distributed `MPIFile` using the given `args` on the given `worker` process.

See `MPIFile`
"""
abstract type DMPIFile{T <: MPIFile} <: MPIFile end
worker(dmdf::DMPIFile) = error("$(typeof(dmdf)) must implement `MPIFiles.worker`")

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
