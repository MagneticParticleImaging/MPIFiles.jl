# Measurements

The low level interface allows to load measured MPI data via the [`measData`](@ref) function.
The returned data is exactly how it is stored on disc. This has the disadvantage
that the user needs to handle the different processing stages in which the data can be stored,
e.g. in time or frequency domain, with or without background correction, or with a permuted frame order.
To cope with this issue, MPIFiles.jl also has a high level interface for loading measurement data.

[`getMeasurements`](@ref) loads the data in time domain. It allows to select frames,
to apply block averaging over frames and periods, to group periods,
and to correct the data for the background and the transfer function of the receive chain.

```julia
f = MPIFile("measurements.mdf")
u = getMeasurements(f, frames=1:100, numAverages=100, bgCorrection=true)
```

[`getMeasurementsFD`](@ref) has the same parameters but returns the data in frequency domain.
Additionally, it can load only a subset of frequencies, e.g. selected with the
[Frequency Filter](@ref), and return the data in a layout suitable for solvers.

```julia
freq = filterFrequencies(f, minFreq=80e3)
uFD = getMeasurementsFD(f, frequencies=freq, numAverages=100)
```

```@docs
getMeasurements
getMeasurementsFD
```
