export ContinuousMechanicalRotationChannel

"Mechanical channel with a continuous rotation."
Base.@kwdef struct ContinuousMechanicalRotationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Frequency of the mechanical rotation. If defined as a vector, the frequencies can swept."
  divider::Union{Integer, Vector{Integer}}
  "Phase of the mechanical rotation. If defined as a vector, the phases will be swept alongside with the frequencies."
  phase::Union{typeof(1.0u"rad"), Vector{typeof(1.0u"rad")}}
end

channeltype(::Type{<:ContinuousMechanicalRotationChannel}) = ContinuousTxChannel()

function createFieldChannel(channelID::AbstractString, channelType::Type{ContinuousMechanicalRotationChannel}, channelDict::Dict{String, Any})
  divider = Int64.(channelDict["divider"])
  phase = uconvert.(u"rad", uparse.(channelDict["phase"]))
  return ContinuousMechanicalRotationChannel(id=channelID, divider=divider, phase=phase)
end