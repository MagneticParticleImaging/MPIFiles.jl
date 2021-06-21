using MPIFiles

filename = "Sequence.toml"
path = joinpath(pwd(), "test", "Sequence", filename)
@info path

sequence = sequenceFromTOML(path)
@info sequence