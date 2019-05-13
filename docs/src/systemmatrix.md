# System Matrices

For loading the system matrix, one could in principle again call `measData` but there
is again a high level function for this job. Since system functions can be very large
it is crutial to load only the subset of frequencies that are used during reconstruction
The high level system matrix loading function is called `getSystemMatrix` and has
the following interface
```julia
function getSystemMatrix(f::MPIFile,
                         frequencies=1:rxNumFrequencies(f)*rxNumChannels(f);
                         bgCorrection=false,
                         loadasreal=false,
                         kargs...)
```
`loadasreal` can again be used when using a solver requiring real numbers.
The most important parameter is frequencies, which defaults to all possible
frequencies over all receive channels. In practice one will determine the
frequencies using the the [Frequency Filter](@ref) functionality.
