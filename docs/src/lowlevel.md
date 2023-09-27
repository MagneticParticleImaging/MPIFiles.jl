# Low Level Interface

The low level interface of MPIFiles.jl consists of a collection of methods that
need to be implemented for each file format. It consists of the following methods
```julia
# general
version, uuid

# study parameters
studyName, studyNumber, studyUuid, studyDescription

# experiment parameters
experimentName, experimentNumber, experimentUuid, experimentDescription,
experimentSubject, experimentIsSimulation, experimentIsCalibration,
experimentHasMeasurement, experimentHasReconstruction

# tracer parameters
tracerName, tracerBatch, tracerVolume, tracerConcentration, tracerSolute,
tracerInjectionTime, tracerVendor

# scanner parameters
scannerFacility, scannerOperator, scannerManufacturer, scannerName, scannerTopology

# acquisition parameters
acqStartTime, acqNumFrames, acqNumAverages, acqGradient, acqOffsetField,
acqNumPeriodsPerFrame

# drive-field parameters
dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform, dfDivider,
dfWaveform, dfCycle

# receiver parameters
rxNumChannels, rxBandwidth, rxNumSamplingPoints, rxTransferFunction, rxUnit,
rxDataConversionFactor, rxInductionFactor

# measurements
measData, measDataTDPeriods, measIsFourierTransformed, measIsTFCorrected,
measIsBGCorrected, measIsFastFrameAxis, measIsFramePermutation, measIsFrequencySelection,
measIsBGFrame, measIsSpectralLeakageCorrected, measFramePermutation

# calibrations
calibSNR, calibFov, calibFovCenter, calibSize, calibOrder, calibPositions,
calibOffsetFields, calibDeltaSampleSize, calibMethod, calibIsMeanderingGrid

# reconstruction results
recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

# additional functions that should be implemented by an MPIFile
filepath, systemMatrixWithBG, systemMatrix, selectedChannels
```
The interface is structured in a similar way as the parameters within the [MDF](https://github.com/MagneticParticleImaging/MDF). Basically, there is a direct mapping between the MDF parameters
and the MPIFiles interface. For instance the parameter `acqNumAvarages` maps to the MDF parameter `/acquisition/numAverages`. Also the dimensionality of the parameters described in the [MDF](https://github.com/MagneticParticleImaging/MDF) is preserved. Thus, the MDF specification can be used as
a documentation of the low level interface of MPIFiles.

## In-memory MDF
Since during data acquisition it is not advisable to write the data directly to an MDF file, there is also an in-memory representation. It can be created either empty by `mdf = MDFv2InMemory()` or with default values set by `mdf = defaultMDFv2InMemory()`. The same interface as described above applies to the in-memory version. When writing it to disk using `saveasMDF("your/file/path/file.mdf", mdf)` it is checked for consistency. If this check should fail with an error, use
`saveasMDF("your/file/path/file.mdf", mdf, failOnInconsistent=true)`. The default is to not fail in order
to not lose measurement data in case of an error.

!!! note
    Note that the dimensions in the MDF documentation are flipped compared to the dimensions in
    Julia. This is because Julia stores the data in column major order, while HDF5 considers row
    major order
