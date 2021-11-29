export StepwiseElectricalChannel

"Electrical channel with a stepwise definition of values."
Base.@kwdef struct StepwiseElectricalChannel <: AcyclicElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Divider of the component."
  divider::Integer
  "values corresponding to the individual steps."
  values::Union{Vector{typeof(1.0u"T")}, Vector{typeof(1.0u"A")}}
end

channeltype(::Type{<:StepwiseElectricalChannel}) = StepwiseTxChannel()

function createFieldChannel(channelID::AbstractString, channelType::Type{StepwiseElectricalChannel}, channelDict::Dict{String, Any})
  divider = channelDict["divider"]
  values = uparse.(channelDict["values"])
  if eltype(values) <: Unitful.Current
    values = values .|> u"A"
  elseif eltype(values) <: Unitful.BField
    values = values .|> u"T"
  else
    error("The values have to be either given as a current or in tesla. You supplied the type `$(eltype(values))`.")
  end

  if mod(divider, length(values)) != 0
    error("The divider $(divider) needs to be a multiple of the $(length(values))")
  end

  return StepwiseElectricalChannel(;id=channelID, divider, values)
end

values(channel::StepwiseElectricalChannel) = channel.values