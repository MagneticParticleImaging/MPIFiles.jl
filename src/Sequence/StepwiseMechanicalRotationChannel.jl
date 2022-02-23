export StepwiseMechanicalRotationChannel

"Mechanical channel with a triggered stepwise rotation."
Base.@kwdef struct StepwiseMechanicalRotationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Step angle of the mechanical rotation."
  stepAngle::typeof(1.0u"rad") #TODO: Should we have a stepsPerSecond value?
end

channeltype(::Type{<:StepwiseMechanicalRotationChannel}) = StepwiseTxChannel()
mechanicalMovementType(::Type{<:StepwiseMechanicalRotationChannel}) = RotationTxChannel()

function createFieldChannel(channelID::AbstractString, channelType::Type{StepwiseMechanicalRotationChannel}, channelDict::Dict{String, Any})
  stepAngle = uconvert.(u"rad", uparse.(channelDict["stepAngle"]))
  return StepwiseMechanicalRotationChannel(id=channelID, stepAngle=stepAngle)
end

cycleDuration(channel::StepwiseMechanicalRotationChannel, baseFrequency::typeof(1.0u"Hz")) = nothing
stepsPerCycle(channel::StepwiseMechanicalRotationChannel) = round(Int64, 2π/channel.stepAngle)