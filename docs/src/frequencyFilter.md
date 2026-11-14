# Frequency Filter

In MPI reconstruction usually only a subset of the frequency components is used, e.g. the
components with a sufficiently large signal-to-noise ratio that are not disturbed by the
excitation. [`filterFrequencies`](@ref) determines such a subset:

```julia
fSM = MPIFile("systemMatrix.mdf")
freq = filterFrequencies(fSM, SNRThresh=5, minFreq=80e3, recChannels=1:2)
```

The SNR threshold requires the SNR estimate of a calibration file, see [`calibSNR`](@ref).
The returned indices can be directly passed to [`getMeasurementsFD`](@ref) and
[`getSystemMatrix`](@ref).

```@docs
filterFrequencies
```
