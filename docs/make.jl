using Documenter, MPIFiles

makedocs(
    format = Documenter.HTML(prettyurls = false),
    modules = [MPIFiles],
    sitename = "MPI Files",
    authors = "Tobias Knopp et al.",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Low Level Interface" => "lowlevel.md",
        "Conversion" => "conversion.md",
        "Measurements" => "measurements.md",
        "System Matrices" => "systemmatrix.md",
        "Frequency Filter" => "frequencyFilter.md",
        "Reconstruction Results" => "reconstruction.md"
      #  "Positions" => "positions.md"
    ],
    warnonly = [:missing_docs]
)

deploydocs(
    repo = "github.com/MagneticParticleImaging/MPIFiles.jl.git",
    target = "build",
)
