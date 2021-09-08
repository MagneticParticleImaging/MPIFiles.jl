@testset "Testing Sequence submodule" begin

filename = "ElectricalSequence.toml"
path = joinpath(@__DIR__, filename)

sequence = sequenceFromTOML(path)

@test sequence.name == "ElectricalSequence"
@test dfDivider(sequence) == [4864 4800; 4864 4800]
@test dfCycle(sequence) ≈ 0.0029184u"s"

end