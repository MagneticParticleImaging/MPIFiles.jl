export StepwiseMechanicalRotationChannel

"Mechanical channel with a triggered stepwise rotation."
Base.@kwdef struct StepwiseMechanicalRotationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Step angle of the mechanical rotation. If defined as a vector, the steps can be non-equidistant."
  stepAngle::Union{typeof(1.0u"rad"), Vector{typeof(1.0u"rad")}}
end

channeltype(::Type{<:StepwiseMechanicalRotationChannel}) = StepwiseTxChannel()

function createFieldChannel(channelID::AbstractString, channelType::Type{StepwiseMechanicalRotationChannel}, channelDict::Dict{String, Any})
  stepAngle = uconvert.(u"rad", uparse.(channelDict["stepAngle"]))
  return StepwiseMechanicalRotationChannel(id=channelID, stepAngle=stepAngle)
end