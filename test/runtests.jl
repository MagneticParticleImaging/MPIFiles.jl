using Dates
using HTTP
using HDF5
using LinearAlgebra
using MPIFiles
using Statistics
using Test
using UUIDs
using Unitful
using Scratch
using LazyArtifacts
using Aqua
using Dagger

const datadir = joinpath(artifact"data", "data")
@info "The test data is located at $datadir."

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."

mkpath(joinpath(tmpdir,"mdf"))
mkpath(joinpath(tmpdir,"mdfim"))
mkpath(joinpath(tmpdir,"positions"))
mkpath(joinpath(tmpdir,"transferFunction"))

@testset "Aqua" begin
  Aqua.test_all(MPIFiles)
end

include("General.jl")
include("Cartesian.jl")
include("Positions.jl")
include("MDFv1.jl")
include("MultiMPIFile.jl")
include("Reco.jl")
include("IMT.jl")
include("CustomSFMeas.jl")
include("MDFInMemory.jl")
include("TransferFunction.jl")
include("FrequencyFilter.jl")
include("MagneticFieldMeasurement.jl")
include("DatasetStore.jl")

@info "The unit tests are done!"
