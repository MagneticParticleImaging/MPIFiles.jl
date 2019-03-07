import ImageMagick
using Documenter, MPIReco

makedocs(
    modules = [MPIReco],
    format = :html,
    sitename = "Julia MPI Package",
    authors = "Tobias Knopp, et al.",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Acquisition Data" => "acquisitionData.md",
        "File Handling" => "filehandling.md",
        "Images" => "image.md",
        "Offresonance" => "offresonance.md",
        "Parallel Imaging" => "SENSE.md",
        "Trajectory" => "trajectories.md",
        "Imaging Operators" => "operators.md",
        "Simulation" => "simulation.md",
        "Reconstruction" => "reconstruction.md",
    ],
    html_prettyurls = false, #!("local" in ARGS),
)

deploydocs(
    repo = "github.com/MagneticParticleImaging/MPIReco.jl.git",
    target = "build",
)
