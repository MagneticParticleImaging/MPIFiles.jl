using Dates
using HTTP
using HDF5
using LinearAlgebra
using MPIFiles
using Statistics
using Test
using UUIDs
using Unitful

if !isdir("data")
  @info "download data.zip"
  HTTP.open("GET", "http://media.tuhh.de/ibi/MPIFiles/data.zip") do http
    open("data.zip", "w") do file
        write(file, http)
    end
  end
  @info "extracting data"
  run(`unzip -oq data.zip`)
  rm("data.zip")
end

include("Positions.jl")
include("General.jl")
include("MDFv1.jl")
include("MultiMPIFile.jl")
include("Reco.jl")
include("IMT.jl")
include("TransferFunction.jl")

@info "The unit tests are done!"
