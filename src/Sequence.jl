using TOML

export Waveform, WAVEFORM_SINE, WAVEFORM_SQUARE, WAVEFORM_TRIANGLE, WAVEFORM_SAWTOOTH_RISING,
       WAVEFORM_SAWTOOTH_FALLING, toWaveform, TxChannel, ElectricalTxChannel, StepwiseElectricalTxChannel,
       MechanicalTxChannel, ElectricalComponent, PeriodicElectricalComponent,
       SweepElectricalComponent, PeriodicElectricalChannel, EquidistantStepwiseElectricalChannel,
       NonequidistantStepwiseElectricalChannel, MechanicalTranslationChannel,
       StepwiseMechanicalRotationChannel, ContinuousMechanicalRotationChannel,
       MagneticField, RxChannel, Sequence, sequenceFromTOML, fieldDictToFields,
       electricalTxChannels, mechanicalTxChannels

"Enum describing the existing waveforms."
@enum Waveform begin
  WAVEFORM_SINE
  WAVEFORM_SQUARE
  WAVEFORM_TRIANGLE
  WAVEFORM_SAWTOOTH_RISING
  WAVEFORM_SAWTOOTH_FALLING
end

waveformRelations = Dict{String, Waveform}(
  "sine" => WAVEFORM_SINE,
  "square" => WAVEFORM_SQUARE,
  "triangle" => WAVEFORM_TRIANGLE,
  "sawtooth_rising" => WAVEFORM_SAWTOOTH_RISING,
  "sawtooth" => WAVEFORM_SAWTOOTH_RISING, # Alias
  "sawtooth_falling" => WAVEFORM_SAWTOOTH_FALLING,
)

toWaveform(value::AbstractString) = waveformRelations[value]

abstract type TxChannel end
abstract type ElectricalTxChannel <: TxChannel end
abstract type StepwiseElectricalTxChannel <: ElectricalTxChannel end
abstract type MechanicalTxChannel <: TxChannel end

abstract type ElectricalComponent end

"Component of an electrical channel with periodic base function."
Base.@kwdef struct PeriodicElectricalComponent <: ElectricalComponent
  "Divider of the component."
  divider::Int64
  "Amplitude (peak) of the component for each period of the field. If defined in Tesla, the calibration configured
  in the scanner will be used."
  amplitude::Vector{Union{typeof(1.0u"T"), typeof(1.0u"V")}} # Is it really the right choice to have the periods here? Or should it be moved to the MagneticField?
  "Phase of the component for each period of the field."
  phase::Vector{typeof(1.0u"rad")}
end

"Sweepable component of an electrical channel with periodic base function."
Base.@kwdef struct SweepElectricalComponent <: ElectricalComponent
  "Divider of the component."
  divider::Vector{Int64}
  "Amplitude (peak) of the channel for each divider in the sweep. If defined in Tesla, the calibration configured
  in the scanner will be used. Must have the same dimension as `dividers`."
  amplitude::Vector{Union{typeof(1.0u"T"), typeof(1.0u"V")}}
  "Phase of the component for each divider in the sweep. If defined as a vector, it must have the same length
  as `dividers`."
  phase::Vector{typeof(1.0u"rad")}
end

"Electrical channel based on based on periodic base functions."
Base.@kwdef struct PeriodicElectricalChannel <: ElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Components added for this channel."
  components::Vector{ElectricalComponent}
  "Waveform of the channel."
  waveform::Waveform = WAVEFORM_SINE
  "Offset of the channel. If defined in Tesla, the calibration configured in the scanner will be used."
  offset::Union{typeof(1.0u"T"), typeof(1.0u"V")} = 0.0u"T"
end

"Electrical channel with a equidistant stepwise definition of DC amplitudes."
Base.@kwdef struct EquidistantStepwiseElectricalChannel <: StepwiseElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Number of steps per field cycle for every equidistant step."
  stepsPerCycle::Int64
  "Amplitudes corresponding to the individual steps."
  amplitude::Vector{Union{typeof(1.0u"T"), typeof(1.0u"V")}}
end

"Electrical channel with a non-equidistant stepwise definition of DC amplitudes."
Base.@kwdef struct NonequidistantStepwiseElectricalChannel <: StepwiseElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Amplitudes corresponding to the individual steps defined in a CSV file."
  sampleFilename::AbstractString
end

"Mechanical channel describing a translational movement."
Base.@kwdef struct MechanicalTranslationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Speed of the channel. If defined as a vector, this must have a length of length(positions)-1."
  speed::Union{typeof(1.0u"m/s"), Vector{typeof(1.0u"m/s")}}
  "Positions that define the endpoints of the movement. The movement is
  repeated after reaching the second endpoint. If the vector contains
  more than two positions, the positions can be moved to in a step-wise
  fashion by using triggers."
  positions::Vector{typeof(1.0u"m")}
end

"Mechanical channel with a stepwise rotation."
Base.@kwdef struct StepwiseMechanicalRotationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Step angle of the mechanical rotation. If defined as a vector, the steps can be non-equidistant."
  stepAngle::Union{typeof(1.0u"rad"), Vector{typeof(1.0u"rad")}}
end

"Mechanical channel with a continuous rotation."
Base.@kwdef struct ContinuousMechanicalRotationChannel <: MechanicalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Frequency of the mechanical rotation. If defined as a vector, the frequencies can swept."
  divider::Union{Int64, Vector{Int64}}
  "Phase of the mechanical rotation. If defined as a vector, the phases will be swept alongside with the frequencies."
  phase::Union{typeof(1.0u"rad"), Vector{typeof(1.0u"rad")}}
end

"""
Description of a magnetic field.

The field can either be electromagnetically or mechanically changed.
The mechanical movement of e.g. an iron yoke would be defined within
two channels, one electrical and one mechanical.
"""
Base.@kwdef struct MagneticField
  "Unique ID of the field description."
  id::AbstractString
  "Transmit channels that are used for the field."
  channels::Vector{TxChannel}

  "Flag if the start of the field should be convoluted. Not used for mechanical fields."
  safeStart::Bool = true
  "Flag if the end of the field should be convoluted. In case of an existing brake on
  a mechanical channel this means a use of the brake."
  safeEnd::Bool = true
  "Flag if the field should be convoluted down in case of an error. In case of an
  existing brake on a mechanical channel this means a use of the brake."
  safeError::Bool = true

  "Flag if the channels of the field should be controlled."
  control::Bool = true
  "Flag if the field should be decoupled. Not used for mechanical channels."
  decouple::Bool = true
end

"Receive channel reference that should be included in the acquisition."
Base.@kwdef struct RxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
end

"""
Description of a sequence that can be run by a scanner.

The sequence can either be continuous or triggered. Triggered in
this context means that the acquisition is done on a certain event,
e.g. the move of a robot. The sweeping of frequencies or movement points
can also be done in a triggered or continuous fashion.
"""
Base.@kwdef struct Sequence
  "Name of the sequence to identify it."
  name::AbstractString
  "Description of the sequence."
  description::AbstractString
  "The scanner targeted by the sequence."
  targetScanner::AbstractString
  "Base frequency for all channels. Mechanical channels are synchronized
  with the electrical ones by referencing the time required for the movement
  against this frequency."
  baseFrequency::typeof(1.0u"Hz")
  "Flag if the sequence has a continuous or triggered acquisition."
  triggered::Bool = false

  "Magnetic fields defined by the sequence."
  fields::Vector{MagneticField}

  "Number of frames to acquire. If `triggered` is true, this number of frames
  is acquired on every trigger. If `triggered` is false, but the dividers
  are defined as a vector in an electrical channel, every frequency will be
  applied for `numFrames`."
  numFrames::Int64 = 1
  "Number of periods within a frame."
  numPeriodsPerFrame::Int64 = 1

  "Receive channels that are used in the sequence."
  rxChannels::Vector{RxChannel}
end

function sequenceFromTOML(filename::AbstractString)
  sequenceDict = TOML.parsefile(filename)

  general = sequenceDict["General"]
  receive = sequenceDict["Receive"]

  splattingDict = Dict{Symbol, Any}()

  splattingDict[:name] = general["name"]
  splattingDict[:description] = general["description"]
  splattingDict[:targetScanner] = general["targetScanner"]
  splattingDict[:baseFrequency] = uparse(general["baseFrequency"])
  if haskey(general, "triggered")
    splattingDict[:triggered] = general["triggered"]
  end
  
  splattingDict[:fields] = fieldDictToFields(sequenceDict["Fields"])

  if haskey(general, "numFrames")
    splattingDict[:numFrames] = receive["numFrames"]
  end
  if haskey(general, "numPeriodsPerFrame")
    splattingDict[:numPeriodsPerFrame] = receive["numPeriodsPerFrame"]
  end

  splattingDict[:rxChannels] = RxChannel.(receive["rxChannels"])

  return Sequence(;splattingDict...)
end

function fieldDictToFields(fieldsDict::Dict{String, Any})
  fields = Vector{MagneticField}()

  rootFields = ["safeStart", "safeEnd", "safeError", "control", "decouple"] # Is reflexion better here?
  for (fieldID, fieldDict) in fieldsDict
    splattingDict = Dict{Symbol, Any}()
    channels = Vector{TxChannel}()
    for (channelID, channelDict) in fieldDict
      # Ignore fields from MagneticField root to get the channels
      if !(channelID in rootFields)
        push!(channels, createFieldChannel(channelID, channelDict))
      else
        splattingDict[Symbol(channelID)] = channelDict
      end
    end
    splattingDict[:id] = fieldID
    splattingDict[:channels] = channels

    field = MagneticField(;splattingDict...)
    push!(fields, field)
  end

  return fields
end

function createFieldChannel(channelID::AbstractString, channelDict::Dict{String, Any})
  # If the channel has further dicts this means it has components and therefore is a PeriodicElectricalChannel
  if any([v isa Dict for (k, v) in channelDict])
    # PeriodicElectricalChannel has optional fields, therefore we go the splatting route
    splattingDict = Dict{Symbol, Any}()
    splattingDict[:id] = channelID

    if haskey(channelDict, "waveform")
      splattingDict[:waveform] = toWaveform(channelDict["waveform"])
    end

    if haskey(channelDict, "offset")
      splattingDict[:offset] = uparse.(channelDict["offset"]) .|> u"T" # TODO: Fails for Volt!
    end

    splattingDict[:components] = Vector{ElectricalComponent}()
    components = [v for (k, v) in channelDict if v isa Dict]
    
    for component in components
      divider = component["divider"]
      amplitude = uparse.(component["amplitude"]) .|> u"T" # TODO: Fails for Volt!
      
      if haskey(component, "phase")
        phase = uparse.(component["phase"])
      else
        phase = repeat(0.0u"rad", length(divider)) # Default phase
      end

      @assert length(amplitude) == length(phase) "The length of amplitude and phase must match."

      if divider isa Vector
        push!(splattingDict[:components],
              SweepElectricalComponent(divider=divider,
                                       amplitude=amplitude,
                                       phase=phase))
      else
        push!(splattingDict[:components],
              PeriodicElectricalComponent(divider=divider,
                                          amplitude=amplitude,
                                          phase=phase))
      end
    end
    return PeriodicElectricalChannel(;splattingDict...)
  elseif haskey(channelDict, "stepsPerCycle") && haskey(channelDict, "amplitude")
    stepsPerCycle = channelDict["stepsPerCycle"]
    amplitude = uparse.(channelDict["amplitude"])

    return EquidistantStepwiseElectricalChannel(id=channelID, stepsPerCycle=stepsPerCycle, amplitude=amplitude)
  elseif haskey(channelDict, "sampleFilename")
    sampleFilename = channelDict["sampleFilename"]

    return NonequidistantStepwiseElectricalChannel(id=channelID, sampleFilename=sampleFilename)
  elseif haskey(channelDict, "speed") && haskey(channelDict, "positions")
    speed = uparse(channelDict["speed"])
    positions = uparse.(channelDict["positions"])

    return MechanicalTranslationChannel(id=channelID, speed=speed, positions=positions)
  elseif haskey(channelDict, "stepAngle")
    stepAngle = uconvert.(u"rad", uparse.(channelDict["stepAngle"]))

    return StepwiseMechanicalRotationChannel(id=channelID, stepAngle=stepAngle)
  elseif haskey(channelDict, "divider") && haskey(channelDict, "phase")
    divider = Int64.(channelDict["divider"])
    phase = uconvert.(u"rad", uparse.(channelDict["phase"]))

    return ContinuousMechanicalRotationChannel(id=channelID, divider=divider, phase=phase)
  else
    @warn "Could not determine type of channel. Channel dict: $channelDict"
  end
end

electricalTxChannels(sequence::Sequence)::Vector{ElectricalTxChannel} = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: ElectricalTxChannel]
mechanicalTxChannels(sequence::Sequence)::Vector{MechanicalTxChannel} = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: MechanicalTxChannel]

