# Function reference

This page collects functions that are not covered by the other pages. The parameters of the
[Low Level Interface](@ref) and the high level functions for loading
[Measurements](@ref), [System Matrices](@ref) and [Reconstruction Results](@ref)
are documented on their respective pages. An index of all documented functions and types is
given at the end of this page.

## File handles
All file handles are subtypes of `MPIFile` and implement the [Low Level Interface](@ref).
Usually, files are opened with the generic [`MPIFile`](@ref) constructor, which selects the
matching type.
```@docs
MPIFile
MDFFile
MDFFileV1
MDFFileV2
BrukerFile
IMTFile
MultiMPIFile
DMPIFile
```

## Derived parameters
These functions are computed from the parameters of the [Low Level Interface](@ref) and are
available for all file formats.
```@docs
acqNumFGFrames
acqNumBGFrames
acqNumPeriods
acqNumPatches
acqNumPeriodsPerPatch
acqFramePeriod
acqFov
acqGradientDiag
acqOffsetFieldShift
rxNumFrequencies
rxFrequencies
rxTimePoints
measFGFrameIdx
measBGFrameIdx
measBGFrameBlockLengths
```

## Mixing factors
```@docs
mixFactorToFreq
mixFactorToFreqIdx
mixingFactors
```

## Positions
```@autodocs
Modules = [MPIFiles]
Pages = ["Positions/Positions.jl", "Positions/Subsampling.jl", "Positions/Interpolation.jl"]
```

## Index
```@index
```
