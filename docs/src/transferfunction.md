# Transfer functions

MPIFiles defines the type `TransferFunction` which represents the system-properties of a linear time-invariant system in frequency space ([see also](https://en.wikipedia.org/wiki/Transfer_function)). In the context of MPIFiles this usually includes the properties of the receive chain of an MPI system, but can also hold information about other parts of the signal chain. 

## Basic construction and data access
The `TransferFunction` object is constructed from samples in frequency space and offers two ways to access the underlying data. The first uses indexing to directly access the underlying data:

```@setup tf
using MPIFiles
using MPIFiles.Unitful
```
```@repl tf
    f = collect(range(0,1e6,step=1e3));
    # TransferFunction of a simple lowpass filter
    tf = TransferFunction(f, 1 ./ (1 .+ im*f/1e4 )) 
    tf[1] # Value of tf at first frequency sample (f = 0 Hz)
```

The second method has to be called on the object itself and uses linear interpolation to access the tf at any frequency:
```@repl tf
    tf(0) # Value of tf at f = 0 Hz
    tf(1e4) # Value of tf at f = 10000 Hz
```

A `TransferFunction` can have multiple channels, which adds a second index to both data access functions. Directly accessing multiple channels is also possible. The complex array used to construct the `TransferFunction` needs to have the shape [number of frequency samples, channels].

```@repl tf
tf = TransferFunction(f, [1 ./(1 .+im*f/1e4) 1 ./(1 .+im*f/1e3)])
tf[11,1]
tf[11,2]
tf(1e4,1)
tf(1e4,2)
tf(1e4, [1,2])
tf(1e4,:)
```
## Units 
To attach units to the `TransferFunction` the keyword-parameter `units` can to be used to give a `Unitful` unit to every channel of the tf. Alternatively `data` can just be a Unitful.Quantity. Then `units` is ignored.

This can be useful if the transfer function is not dimensionless but relates two physical quantities, e.g. voltage and current in the form of an impedance. All **interpolated** accesses to tf data then return a `Unitful.Quantity`.

```@repl tf
R = 1; # Ohm
L = 10e-6; # Henry
RL = TransferFunction(f, R .+ im*2pi*f*L, units=["V/A"])
RL([0,100e3])
```

```@repl tf
f_unitful = collect(range(0u"Hz",1u"MHz",step=1u"kHz"));
R = 1u"Ω";
L = 10u"µH";
RL = TransferFunction(f_unitful, R .+ im*2pi*f_unitful*L .|> u"V/A")
RL([0,100e3])
```

## Saving and loading
A `TransferFunction` object can be saved to and loaded from a .h5 file.

```@docs
MPIFiles.save(filename::String, tf::TransferFunction)
MPIFiles.TransferFunction(::String)
```

## Constructors
The `TransferFunction` constructor can take either a complex data array or two arrays representing the amplitude and phase of the transfer function. Unitful conversion is automatically done for all parameters.

It is also possible to construct a `TransferFunction` from the transfer function data included in an `MPIFile`.

```@docs
    TransferFunction
    MPIFiles.TransferFunction(::MPIFile)
```

## Other interesting functions
```@docs
  combine(::TransferFunction, ::TransferFunction)
  MPIFiles.load_tf_fromVNA
  MPIFiles.processRxTransferFunction
```

