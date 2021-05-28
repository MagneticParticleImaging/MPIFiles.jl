using TOML

export Waveform, WAVEFORM_SINE, WAVEFORM_SQUARE, WAVEFORM_TRIANGLE, WAVEFORM_SAWTOOTH_RISING,
       WAVEFORM_SAWTOOTH_FALLING, toWaveform, fromWaveform, TxChannel, ElectricalTxChannel, StepwiseElectricalTxChannel,
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
fromWaveform(value::Waveform) = [k for (k, v) in waveformRelations if v == value][1]

abstract type TxChannel end
abstract type ElectricalTxChannel <: TxChannel end
abstract type StepwiseElectricalTxChannel <: ElectricalTxChannel end
abstract type MechanicalTxChannel <: TxChannel end

abstract type ElectricalComponent end

"Component of an electrical channel with periodic base function."
Base.@kwdef struct PeriodicElectricalComponent <: ElectricalComponent
  "Divider of the component."
  divider::Integer
  "Amplitude (peak) of the component for each period of the field. If defined in Tesla, the calibration configured
  in the scanner will be used."
  amplitude::Vector{Union{typeof(1.0u"T"), typeof(1.0u"V")}} # Is it really the right choice to have the periods here? Or should it be moved to the MagneticField?
  "Phase of the component for each period of the field."
  phase::Vector{typeof(1.0u"rad")}
  "Waveform of the component."
  waveform::Waveform = WAVEFORM_SINE
end

"Sweepable component of an electrical channel with periodic base function."
Base.@kwdef struct SweepElectricalComponent <: ElectricalComponent
  "Divider of the component."
  divider::Vector{Integer}
  "Amplitude (peak) of the channel for each divider in the sweep. If defined in Tesla, the calibration configured
  in the scanner will be used. Must have the same dimension as `dividers`."
  amplitude::Vector{Union{typeof(1.0u"T"), typeof(1.0u"V")}}
  "Phase of the component for each divider in the sweep. If defined as a vector, it must have the same length
  as `dividers`."
  phase::Vector{typeof(1.0u"rad")}
  "Waveform of the component."
  waveform::Waveform = WAVEFORM_SINE
end

"Electrical channel based on based on periodic base functions."
Base.@kwdef struct PeriodicElectricalChannel <: ElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Components added for this channel."
  components::Vector{ElectricalComponent}
  "Offset of the channel. If defined in Tesla, the calibration configured in the scanner will be used."
  offset::Union{typeof(1.0u"T"), typeof(1.0u"V")} = 0.0u"T"
end

"Electrical channel with a equidistant stepwise definition of DC amplitudes."
Base.@kwdef struct EquidistantStepwiseElectricalChannel <: StepwiseElectricalTxChannel
  "ID corresponding to the channel configured in the scanner."
  id::AbstractString
  "Number of steps per field cycle for every equidistant step."
  stepsPerCycle::Integer
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

"Mechanical channel with a triggered stepwise rotation."
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
  divider::Union{Integer, Vector{Integer}}
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

  "Flag if the start of the field should be convoluted.
  If the DAQ does not support this, it can may fall back
  to postponing the application of the settings.
  Not used for mechanical fields."
  safeStart::Bool = true
  "Flag if a transition of the field should be convoluted.
  If the DAQ does not support this, it can may fall back
  to postponing the application of the settings.
  Not used for mechanical fields."
  safeTransition::Bool = true
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
  against this frequency. Please note that the latter has uncertainties."
  txBaseFrequency::typeof(1.0u"Hz")
  "Flag if the sequence has a continuous or triggered acquisition."
  triggered::Bool = false

  "Magnetic fields defined by the sequence."
  fields::Vector{MagneticField}

  "Number of frames to acquire. If `triggered` is true, this number of frames
  is acquired on every trigger. If `numFrames` is a Vector, each trigger will
  acquire the given amount of frames.  If `triggered` is false, but the dividers
  are defined as a vector in an electrical channel, every frequency will be
  applied for `numFrames`."
  numFrames::Union{Integer, Vector{<:Integer}} = 1
  "Number of periods within a frame."
  numPeriodsPerFrame::Integer = 1
  "Number of block averages per period."
  numAverages::Integer = 1
  "Bandwidth (half the sample rate) of the receiver. In DAQs which decimate the data,
  this also determines the decimation. Note: this is currently a
  scalar since the MDF does not allow for multiple sampling rates yet."
  rxBandwidth::typeof(1.0u"Hz")

  "Receive channels that are used in the sequence."
  rxChannels::Vector{RxChannel}
end

function Sequence(filename::AbstractString)
  return sequenceFromTOML(filename)
end

function sequenceFromTOML(filename::AbstractString)
  sequenceDict = TOML.parsefile(filename)

  general = sequenceDict["General"]
  receive = sequenceDict["Receive"]

  splattingDict = Dict{Symbol, Any}()

  splattingDict[:name] = general["name"]
  splattingDict[:description] = general["description"]
  splattingDict[:targetScanner] = general["targetScanner"]
  splattingDict[:txBaseFrequency] = uparse(general["txBaseFrequency"])
  splattingDict[:rxBandwidth] = uparse(general["rxBandwidth"])
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

  sequence =  Sequence(;splattingDict...)

  # TODO: Sanity check on sequence (equal length of triggered vectors etc.)

  return sequence
end

function fieldDictToFields(fieldsDict::Dict{String, Any})
  fields = Vector{MagneticField}()

  rootFields = ["safeStart", "safeTransition", "safeEnd", "safeError", "control", "decouple"] # Is reflexion better here?
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

      if haskey(component, "waveform")
        waveform = toWaveform(component["waveform"])
      else
        waveform = WAVEFORM_SINE # Default to sine
      end

      @assert length(amplitude) == length(phase) "The length of amplitude and phase must match."

      if divider isa Vector
        push!(splattingDict[:components],
              SweepElectricalComponent(divider=divider,
                                       amplitude=amplitude,
                                       phase=phase,
                                       waveform=waveform))
      else
        push!(splattingDict[:components],
              PeriodicElectricalComponent(divider=divider,
                                          amplitude=amplitude,
                                          phase=phase,
                                          waveform=waveform))
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

acqGradient(sequence::Sequence) = nothing # TODO: Implement
acqNumAverages(sequence::Sequence) = sequence.numAverages
acqNumFrames(sequence::Sequence) = sum(sequence.numFrames) # TODO: What about triggered sequences?
acqNumPeriodsPerFrame(sequence::Sequence) = sequence.numPeriodsPerFrame
acqOffsetField(sequence::Sequence) = nothing # TODO: Implement

dfBaseFrequency(sequence::Sequence) = sequence.txBaseFrequency
txBaseFrequency(sequence::Sequence) = sequence.txBaseFrequency
dfCycle(sequence::Sequence) = lcm(dfDivider(sequence))/dfBaseFrequency(sequence)

function dfDivider(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: PeriodicElectricalChannel]
  maxComponents = maximum([length(channel.components) for channel in channels])
  result = zeros(Int64, (dfNumChannels(sequence), maxComponents))
  for (channelIdx, channel) in enumerate(channels)
    for (componentIdx, component) in enumerate(channel.components)
      result[channelIdx, componentIdx] = component.divider
    end
  end
  return result
end

dfNumChannels(sequence::Sequence) = length(electricalTxChannels(sequence)) # TODO: How do we integrate the mechanical channels?

function dfPhase(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: PeriodicElectricalChannel]
  maxComponents = maximum([length(channel.components) for channel in channels])
  numPeriods = length(channels[1].components[1].phase) # Should all be of the same length
  result = zeros(typeof(1.0u"rad"), (numPeriods, dfNumChannels(sequence), maxComponents))
  for (channelIdx, channel) in enumerate(channels)
    for (componentIdx, component) in enumerate(channel.components)
      for (periodIdx, phase) in enumerate(component.phase)
        result[periodIdx, channelIdx, componentIdx] = phase
      end
    end
  end
  return result
end

function dfStrength(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: PeriodicElectricalChannel]
  maxComponents = maximum([length(channel.components) for channel in channels])
  numPeriods = length(channels[1].components[1].amplitude) # Should all be of the same length
  result = zeros(typeof(1.0u"T"), (numPeriods, dfNumChannels(sequence), maxComponents))
  for (channelIdx, channel) in enumerate(channels)
    for (componentIdx, component) in enumerate(channel.components)
      for (periodIdx, strength) in enumerate(component.amplitude) # TODO: What do we do if this is in volt? The conversion factor is with the scanner... Remove the volt version?
        result[periodIdx, channelIdx, componentIdx] = strength
      end
    end
  end
  return result
end

function dfWaveform(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: PeriodicElectricalChannel]
  maxComponents = maximum([length(channel.components) for channel in channels])
  result = fill(WAVEFORM_SINE, (dfNumChannels(sequence), maxComponents))
  for (channelIdx, channel) in enumerate(channels)
    for (componentIdx, component) in enumerate(channel.components)
      result[channelIdx, componentIdx] = component.waveform
    end
  end
  return result
end

rxBandwidth(sequence::Sequence) = sequence.rxBandwidth
rxNumChannels(sequence::Sequence) = length(sequence.rxChannels)
rxNumSamplingPoints(sequence::Sequence) = round(Int64, ustrip(NoUnits, rxBandwidth(sequence)*2*dfCycle(sequence)))