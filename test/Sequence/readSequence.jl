using MPIFiles

filename = "Sequence.toml"
path = joinpath(pwd(), "scratch", filename)
@info path

sequence = sequenceFromTOML(path)
@info sequence