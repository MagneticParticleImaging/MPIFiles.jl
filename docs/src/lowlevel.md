# Low Level Interface

The low level interface of MPIFiles.jl consists of a collection of methods that
need to be implemented for each file format. The interface is structured in a similar way as the parameters within the [MDF](https://github.com/MagneticParticleImaging/MDF). Basically, there is a direct mapping between the MDF parameters
and the MPIFiles interface. For instance the parameter `acqNumAverages` maps to the MDF parameter `/acquisition/numAverages`. The MDF path of each parameter is given in the
documentation of the corresponding function below.

Parameters that are optional in the MDF return `nothing` if they are not available. Parameters that are required but missing in the file return `missing`.

!!! note
    The dimensions in the MDF specification are flipped compared to the dimensions in
    Julia. This is because Julia stores the data in column major order, while HDF5 considers row
    major order. The sizes given below are in Julia order.

The sizes use the following symbols of the MDF specification:

| Symbol | Meaning |
|:-------|:--------|
| A | number of tracers |
| N | number of acquired frames (N = O + E) |
| O, E | number of foreground and background frames |
| J | number of periods within one frame |
| Y | number of partitions of each patch position |
| C | number of receive channels |
| D | number of drive-field channels |
| F | number of frequencies describing the drive-field waveform |
| V | number of sampling points during one drive-field period |
| K | number of frequencies describing the processed data |
| P, Q, S | number of voxels, frames and channels of the reconstruction |

## General parameters
```@docs
version(::MPIFile)
uuid(::MPIFile)
MPIFiles.time(::MPIFile)
```

## Study parameters
```@docs
studyName(::MPIFile)
studyNumber(::MPIFile)
studyUuid(::MPIFile)
studyDescription(::MPIFile)
studyTime(::MPIFile)
```

## Experiment parameters
```@docs
experimentName(::MPIFile)
experimentNumber(::MPIFile)
experimentUuid(::MPIFile)
experimentDescription(::MPIFile)
experimentSubject(::MPIFile)
experimentIsSimulation(::MPIFile)
experimentIsCalibration(::MPIFile)
experimentHasMeasurement(::MPIFile)
experimentHasReconstruction(::MPIFile)
```

## Tracer parameters
```@docs
tracerName(::MPIFile)
tracerBatch(::MPIFile)
tracerVolume(::MPIFile)
tracerConcentration(::MPIFile)
tracerSolute(::MPIFile)
tracerInjectionTime(::MPIFile)
```

## Scanner parameters
```@docs
scannerBoreSize(::MPIFile)
scannerFacility(::MPIFile)
scannerOperator(::MPIFile)
scannerManufacturer(::MPIFile)
scannerName(::MPIFile)
scannerTopology(::MPIFile)
```

## Acquisition parameters
```@docs
acqStartTime(::MPIFile)
acqNumFrames(::MPIFile)
acqNumAverages(::MPIFile)
acqNumPeriodsPerFrame(::MPIFile)
acqGradient(::MPIFile)
acqOffsetField(::MPIFile)
```

## Drive-field parameters
```@docs
dfNumChannels(::MPIFile)
dfStrength(::MPIFile)
dfPhase(::MPIFile)
dfBaseFrequency(::MPIFile)
dfDivider(::MPIFile)
dfWaveform(::MPIFile)
dfCustomWaveform(::MPIFile)
dfCycle(::MPIFile)
```

## Receiver parameters
```@docs
rxNumChannels(::MPIFile)
rxBandwidth(::MPIFile)
rxNumSamplingPoints(::MPIFile)
rxTransferFunction(::MPIFile)
rxTransferFunctionFileName(::MPIFile)
rxHasTransferFunction(::MPIFile)
rxUnit(::MPIFile)
rxDataConversionFactor(::MPIFile)
rxInductionFactor(::MPIFile)
```

## Measurement parameters
```@docs
measData(::MPIFile)
MPIFiles.measDataTD(::MPIFile)
measDataTDPeriods(::MPIFile, ::Any)
measIsFourierTransformed(::MPIFile)
measIsTFCorrected(::MPIFile)
measIsBGCorrected(::MPIFile)
measIsFastFrameAxis(::MPIFile)
measIsFramePermutation(::MPIFile)
measIsFrequencySelection(::MPIFile)
measIsBGFrame(::MPIFile)
measIsSpectralLeakageCorrected(::MPIFile)
measIsSparsityTransformed(::MPIFile)
measIsCalibProcessed(::MPIFile)
measFramePermutation(::MPIFile)
measFrequencySelection(::MPIFile)
measSparsityTransformation(::MPIFile)
measSubsamplingIndices(::MPIFile)
```

## Calibration parameters
```@docs
calibSNR(::MPIFile)
calibFov(::MPIFile)
calibFovCenter(::MPIFile)
calibSize(::MPIFile)
calibOrder(::MPIFile)
calibPositions(::MPIFile)
calibOffsetFields(::MPIFile)
calibDeltaSampleSize(::MPIFile)
calibMethod(::MPIFile)
calibIsMeanderingGrid(::MPIFile)
```

## Reconstruction parameters
```@docs
recoData(::MPIFile)
recoFov(::MPIFile)
recoFovCenter(::MPIFile)
recoSize(::MPIFile)
recoOrder(::MPIFile)
recoPositions(::MPIFile)
MPIFiles.recoIsOverscanRegion(::MPIFile)
```

## Additional functions
```@docs
filepath(::MPIFile)
```

## In-memory MDF
Since during data acquisition it is not advisable to write the data directly to an MDF file, there is also an in-memory representation. It can be created either empty by `mdf = MDFv2InMemory()` or with default values set by `mdf = defaultMDFv2InMemory()`. The same interface as described above applies to the in-memory version. Additionally, every parameter can be set by calling the accessor function with the new value as second argument, e.g. `acqNumFrames(mdf, 100)`.

When writing it to disk using `saveasMDF("your/file/path/file.mdf", mdf)` it is checked for consistency. If this check should fail with an error, use
`saveasMDF("your/file/path/file.mdf", mdf, failOnInconsistent=true)`. The default is to not fail in order
to not lose measurement data in case of an error.

The in-memory MDF consists of the following parts, which correspond to the groups of the MDF.

```@docs
MDFv2InMemory
defaultMDFv2InMemory
MDFv2Root
MDFv2Study
MDFv2Experiment
MDFv2Tracer
MDFv2Scanner
MDFv2Acquisition
MDFv2Drivefield
MDFv2Receiver
MDFv2Measurement
MDFv2Calibration
MDFv2Reconstruction
MDFv2Variables
```
