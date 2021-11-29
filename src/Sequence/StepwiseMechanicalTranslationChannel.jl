export StepwiseMechanicalTranslationChannel

"Mechanical channel describing a stepwise translational movement."
Base.@kwdef struct StepwiseMechanicalTranslationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Speed of the channel. If defined as a vector, this must have a length of length(positions)-1."
  speed::Union{typeof(1.0u"m/s"), Vector{typeof(1.0u"m/s")}}
  "Positions that define the steps of the movement."
  positions::Vector{typeof(1.0u"m")}
end

channeltype(::Type{<:StepwiseMechanicalTranslationChannel}) = StepwiseTxChannel()

function createFieldChannel(channelID::AbstractString, channelType::Type{StepwiseMechanicalTranslationChannel}, channelDict::Dict{String, Any})
  speed = uparse(channelDict["speed"])
  positions = uparse.(channelDict["positions"])
  return StepwiseMechanicalTranslationChannel(id=channelID, speed=speed, positions=positions)
end