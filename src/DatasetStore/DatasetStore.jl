export DatasetStore, studydir, path, readonly, remove

include("Utils.jl")

abstract type DatasetStore end

@mustimplement studydir(::DatasetStore)
@mustimplement readonly(::DatasetStore)

include("Study.jl")
include("Experiment.jl")
include("MDFDatasetStore.jl")
include("Reconstruction.jl")
include("Visualization.jl")
include("BrukerDatasetStore.jl")
include("Export.jl")
include("SFDatabase.jl")

