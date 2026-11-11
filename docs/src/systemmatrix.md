# System Matrices

For loading the system matrix, one could in principle again call [`measData`](@ref), but there
is again a high level function for this job. Since system matrices can be very large
it is crucial to load only the subset of frequencies that are used during reconstruction.
In practice, one will determine the frequencies using the [Frequency Filter](@ref) functionality:

```julia
fSM = MPIFile("systemMatrix.mdf")
freq = filterFrequencies(fSM, SNRThresh=5, minFreq=80e3)
S = getSystemMatrix(fSM, freq, bgCorrection=true)
```

```@docs
getSystemMatrix
```

The grid on which the system matrix was measured can be obtained with
```@docs
MPIFiles.calibAxis
```
