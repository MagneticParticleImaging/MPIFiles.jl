export DatasetStore, studydir, path, readonly, remove

include("Utils.jl")

abstract type DatasetStore end

# mandatory functions
@mustimplement studydir(::DatasetStore)
@mustimplement readonly(::DatasetStore)

# functions for writable store
@mustimplement Base.empty!(d::DatasetStore)
@mustimplement changeParam(exp, paramName::AbstractString, paramValue)


include("Study.jl")
include("Experiment.jl")
include("MDFDatasetStore.jl")
include("Reconstruction.jl")
include("Visualization.jl")
include("BrukerDatasetStore.jl")
include("Export.jl")
include("SFDatabase.jl")

