using TOML

export Waveform, WAVEFORM_SINE, WAVEFORM_SQUARE, WAVEFORM_TRIANGLE, WAVEFORM_SAWTOOTH_RISING,
       WAVEFORM_SAWTOOTH_FALLING, toWaveform, fromWaveform, TxChannel, ElectricalTxChannel, StepwiseElectricalTxChannel,
       MechanicalTxChannel, ElectricalComponent, PeriodicElectricalComponent,
       SweepElectricalComponent, PeriodicElectricalChannel, EquidistantStepwiseElectricalChannel,
       NonequidistantStepwiseElectricalChannel, MechanicalTranslationChannel,
       StepwiseMechanicalRotationChannel, ContinuousMechanicalRotationChannel,
       MagneticField, RxChannel, AcquisitionSettings, Sequence, sequenceFromTOML, fieldDictToFields,
       id, offset, components, divider, amplitude, phase, waveform, electricalTxChannels,
       mechanicalTxChannels, periodicElectricalTxChannels, stepwiseElectricalTxChannels,
       acqGradient, acqNumFrames, acqNumPeriodsPerFrame,
       acqNumAverages, acqNumFrameAverages, acqOffsetField, numForegroundTriggers, numBackgroundTriggers,
       numTriggers, numTriggersTotal, dfBaseFrequency, txBaseFrequency,
       txCycle, dfDivider, dfNumChannels, dfPhase, dfStrength, dfWaveform, rxBandwidth,
       rxNumChannels, rxNumSamplingPoints, rxNumSamplesPerPeriod, rxChannels,
       needsControl, needsDecoupling, needsControlOrDecoupling

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
  "Amplitude (peak) of the component for each period of the field."
  amplitude::Union{Vector{typeof(1.0u"T")}, Vector{typeof(1.0u"V")}} # Is it really the right choice to have the periods here? Or should it be moved to the MagneticField?
  "Phase of the component for each period of the field."
  phase::Vector{typeof(1.0u"rad")}
  "Waveform of the component."
  waveform::Waveform = WAVEFORM_SINE
end

"Sweepable component of an electrical channel with periodic base function.
Note: Does not allow for changes in phase since this would make the switch
between frequencies difficult."
Base.@kwdef struct SweepElectricalComponent <: ElectricalComponent
  "Divider of the component."
  divider::Vector{Integer}
  "Amplitude (peak) of the channel for each divider in the sweep. Must have the same dimension as `divider`."
  amplitude::Union{Vector{typeof(1.0u"T")}, Vector{typeof(1.0u"V")}}
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
  amplitude::Union{Vector{typeof(1.0u"T")}, Vector{typeof(1.0u"V")}}
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
  safeStartInterval::typeof(1.0u"s") = 0.5u"s"
  "Flag if a transition of the field should be convoluted.
  If the DAQ does not support this, it can may fall back
  to postponing the application of the settings.
  Not used for mechanical fields."
  safeTransitionInterval::typeof(1.0u"s") = 0.5u"s"
  "Flag if the end of the field should be convoluted. In case of an existing brake on
  a mechanical channel this means a use of the brake."
  safeEndInterval::typeof(1.0u"s") = 0.5u"s"
  "Flag if the field should be convoluted down in case of an error. In case of an
  existing brake on a mechanical channel this means a use of the brake."
  safeErrorInterval::typeof(1.0u"s") = 0.5u"s"

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

"Settings for acquisition on a given trigger."
struct AcquisitionTriggerSettings{N}
  "Number of frames to acquire. If `triggered` is true, this number of frames
  is acquired on every trigger. If `numFrames` is a tuple, each trigger will
  acquire the given amount of frames.  If `triggered` is false, but the dividers
  are defined as a vector in an electrical channel, every frequency will be
  applied for `numFrames`."
  numFrames::Union{Integer, NTuple{N, <:Integer}}
  "Number of block averages per period."
  numAverages::Union{Integer, NTuple{N, <:Integer}}
  "Number of frames to average blockwise."
  numFrameAverages::Union{Integer, NTuple{N, <:Integer}}

  # I am not completely sure, why I can't do this directly, but having default
  # values on the fields makes the type definition fail
  function AcquisitionTriggerSettings(;
      numFrames::Union{Integer, NTuple{N, <:Integer}} = 1,
      numAverages::Union{Integer, NTuple{N, <:Integer}} = 1,
      numFrameAverages::Union{Integer, NTuple{N, <:Integer}} = 1
    ) where {N}
    
    N2 = maximum([length(numFrames), length(numAverages), length(numFrameAverages)])
    new{N2}(numFrames, numAverages, numFrameAverages)
  end
end

"Settings for acquiring the sequence."
Base.@kwdef struct AcquisitionSettings
  "Receive channels that are used in the sequence."
  channels::Vector{RxChannel}
  "Bandwidth (half the sample rate) of the receiver. In DAQs which decimate the data,
  this also determines the decimation. Note: this is currently a
  scalar since the MDF does not allow for multiple sampling rates yet."
  bandwidth::typeof(1.0u"Hz")
  "Number of periods within a frame."
  numPeriodsPerFrame::Integer = 1
  
  "Setting for foreground triggers."
  foreground::AcquisitionTriggerSettings = AcquisitionTriggerSettings()
  "Setting for background triggers."
  background::AcquisitionTriggerSettings = AcquisitionTriggerSettings()
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
  baseFrequency::typeof(1.0u"Hz")
  "Flag if the sequence has a continuous or triggered acquisition."
  triggered::Bool = false

  "Magnetic fields defined by the sequence."
  fields::Vector{MagneticField}

  "Settings for the acquisition."
  acquisiton::AcquisitionSettings
end

function Sequence(filename::AbstractString)
  return sequenceFromTOML(filename)
end

function sequenceFromTOML(filename::AbstractString)
  sequenceDict = TOML.parsefile(filename)

  general = sequenceDict["General"]
  acquisition = sequenceDict["Acquisition"]

  splattingDict = Dict{Symbol, Any}()

  # Main section
  splattingDict[:name] = general["name"]
  splattingDict[:description] = general["description"]
  splattingDict[:targetScanner] = general["targetScanner"]
  splattingDict[:baseFrequency] = uparse(general["baseFrequency"])
  if haskey(general, "triggered")
    splattingDict[:triggered] = general["triggered"]
  end
  
  # Fields
  splattingDict[:fields] = fieldDictToFields(sequenceDict["Fields"])

  # Acquisition
  acqSplattingDict = Dict{Symbol, Any}()
  
  acqSplattingDict[:channels] = RxChannel.(acquisition["channels"])
  acqSplattingDict[:bandwidth] = uparse(acquisition["bandwidth"])
  if haskey(acquisition, "numPeriodsPerFrame")
    acqSplattingDict[:numPeriodsPerFrame] = acquisition["numPeriodsPerFrame"]
  end

  # Acquisition: Triggers
  for (symbol_, key) in [(:foreground, "Foreground"), (:background, "Background")]
    triggerPart = acquisition[key]
    acqTriggerSplattingDict = Dict{Symbol, Any}()
    
    if haskey(triggerPart, "numFrames")
      acqTriggerSplattingDict[:numFrames] = triggerPart["numFrames"]
    end
    if haskey(triggerPart, "numAverages")
      acqTriggerSplattingDict[:numAverages] = triggerPart["numAverages"]
    end
    if haskey(triggerPart, "numFrameAverages")
      acqTriggerSplattingDict[:numFrameAverages] = triggerPart["numFrameAverages"]
    end

    acqSplattingDict[symbol_] = AcquisitionTriggerSettings(;acqTriggerSplattingDict...)
  end

  splattingDict[:acquisiton] = AcquisitionSettings(;acqSplattingDict...)

  sequence =  Sequence(;splattingDict...)

  # TODO: Sanity check on sequence (equal length of triggered vectors etc.)

  return sequence
end

function fieldDictToFields(fieldsDict::Dict{String, Any})
  fields = Vector{MagneticField}()

  rootFields = ["safeStartInterval", "safeTransitionInterval", "safeEndInterval", "safeErrorInterval", "control", "decouple"] # Is reflexion better here?
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
    
    if haskey(fieldDict, "safeStartInterval")
      splattingDict[:safeStartInterval] = uparse(fieldDict["safeStartInterval"])
    end
    if haskey(fieldDict, "safeTransitionInterval")
      splattingDict[:safeTransitionInterval] = uparse(fieldDict["safeTransitionInterval"])
    end
    if haskey(fieldDict, "safeEndInterval")
      splattingDict[:safeEndInterval] = uparse(fieldDict["safeEndInterval"])
    end
    if haskey(fieldDict, "safeErrorInterval")
      splattingDict[:safeErrorInterval] = uparse(fieldDict["safeErrorInterval"])
    end

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
      tmp = uparse.(channelDict["offset"])
      if eltype(tmp) <: Unitful.Voltage
        tmp = tmp .|> u"V"
      elseif eltype(tmp) <: Unitful.BField
        tmp = tmp .|> u"T"
      else
        error("The value for an offset has to be either given as a voltage or in tesla. You supplied the type `$(eltype(tmp))`.")
      end
      splattingDict[:offset] = tmp
    end

    splattingDict[:components] = Vector{ElectricalComponent}()
    components = [v for (k, v) in channelDict if v isa Dict]
    
    for component in components
      divider = component["divider"]

      amplitude = uparse.(component["amplitude"])
      if eltype(amplitude) <: Unitful.Voltage
        amplitude = amplitude .|> u"V"
      elseif eltype(amplitude) <: Unitful.BField
        amplitude = amplitude .|> u"T"
      else
        error("The value for an amplitude has to be either given as a voltage or in tesla. You supplied the type `$(eltype(tmp))`.")
      end
      
      if haskey(component, "phase")
        phase = uparse.(component["phase"])
      else
        phase = fill(0.0u"rad", length(divider)) # Default phase
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
    if eltype(amplitude) <: Unitful.Voltage
      amplitude = amplitude .|> u"V"
    elseif eltype(amplitude) <: Unitful.BField
      amplitude = amplitude .|> u"T"
    else
      error("The value for an amplitude has to be either given as a voltage or in tesla. You supplied the type `$(eltype(tmp))`.")
    end

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

name(sequence::Sequence) = sequence.name
description(sequence::Sequence) = sequence.description
targetScanner(sequence::Sequence) = sequence.targetScanner
baseFrequency(sequence::Sequence) = sequence.baseFrequency
isTriggered(sequence::Sequence) = sequence.triggered

electricalTxChannels(sequence::Sequence)::Vector{ElectricalTxChannel} = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: ElectricalTxChannel]
mechanicalTxChannels(sequence::Sequence)::Vector{MechanicalTxChannel} = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: MechanicalTxChannel]
periodicElectricalTxChannels(sequence::Sequence)::Vector{PeriodicElectricalChannel} = [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: PeriodicElectricalChannel]
stepwiseElectricalTxChannels(sequence::Sequence)::Vector{EquidistantStepwiseElectricalChannel} = 
  [channel for field in sequence.fields for channel in field.channels if typeof(channel) <: EquidistantStepwiseElectricalChannel]


id(channel::TxChannel) = channel.id
offset(channel::PeriodicElectricalChannel) = channel.offset

# Periodic components
components(channel::PeriodicElectricalChannel) = channel.components
divider(component::ElectricalComponent, trigger::Integer=1) = length(component.divider) == 1 ? component.divider[1] : component.divider[trigger]
amplitude(component::PeriodicElectricalComponent; period::Integer=1) = component.amplitude[period]
amplitude(component::SweepElectricalComponent; trigger::Integer=1) = component.amplitude[period]
phase(component::PeriodicElectricalComponent, trigger::Integer=1) = component.phase[trigger]
phase(component::SweepElectricalComponent, trigger::Integer=1) = 0.0u"rad"
waveform(component::ElectricalComponent) = component.waveform

acqGradient(sequence::Sequence) = nothing # TODO: Implement
function acqNumPeriodsPerFrame(sequence::Sequence)
  channels = stepwiseElectricalTxChannels(sequence)
  numPeriods = [ c.stepsPerCycle*length(c.amplitude) for c in channels ]
  if minimum(numPeriods) != maximum(numPeriods)
    error("Sequence contains stepwise electrical channels of different length")
  end
  return first(numPeriods)
end
acqForegroundNumFrames(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.foreground.numFrames[trigger] : sequence.acquisiton.foreground.numFrames
acqForegroundNumAverages(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.foreground.numAverages[trigger] : sequence.acquisiton.foreground.numAverages
acqForegroundNumFrameAverages(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.foreground.numFrameAverages[trigger] : sequence.acquisiton.foreground.numFrameAverages
acqBackgroundNumFrames(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.background.numFrames[trigger] : sequence.acquisiton.background.numFrames
acqBackgroundNumAverages(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.background.numAverages[trigger] : sequence.acquisiton.background.numAverages
acqBackgroundNumFrameAverages(sequence::Sequence, trigger::Integer=1) = trigger>1 ? sequence.acquisiton.background.numFrameAverages[trigger] : sequence.acquisiton.background.numFrameAverages
acqNumFrames(sequence::Sequence, trigger::Integer=1, isBackground::Bool=false) = isBackground ? acqBackgroundNumFrames(sequence, trigger) : acqForegroundNumFrames(sequence, trigger)
acqNumAverages(sequence::Sequence, trigger::Integer=1, isBackground::Bool=false) = isBackground ? acqBackgroundNumAverages(sequence, trigger) : acqForegroundNumAverages(sequence, trigger)
acqNumFrameAverages(sequence::Sequence, trigger::Integer=1, isBackground::Bool=false) = isBackground ? acqBackgroundNumFrameAverages(sequence, trigger) : acqForegroundNumFrameAverages(sequence, trigger)
acqOffsetField(sequence::Sequence) = nothing # TODO: Implement

function numForegroundTriggers(sequence::Sequence)
  part = sequence.acquisiton.foreground
  return maximum([length(part.numFrames), length(part.numAverages), length(part.numFrameAverages)])
end
function numBackgroundTriggers(sequence::Sequence)
  part = sequence.acquisiton.background
  return maximum([length(part.numFrames), length(part.numAverages), length(part.numFrameAverages)])
end
numTriggers(sequence::Sequence, isBackground::Bool=false) = isBackground ? numBackgroundTriggers(sequence) : numForegroundTriggers(sequence)
numTriggersTotal(sequence::Sequence) = numForegroundTriggers(sequence)+numBackgroundTriggers(sequence)

dfBaseFrequency(sequence::Sequence) = sequence.baseFrequency
txBaseFrequency(sequence::Sequence) = dfBaseFrequency(sequence) # Alias, since this might not only concern the drivefield
dfCycle(sequence::Sequence) = lcm(dfDivider(sequence))/dfBaseFrequency(sequence) |> u"s"
txCycle(sequence::Sequence) = dfCycle(sequence) # Alias, since this might not only concern the drivefield

function dfDivider(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = periodicElectricalTxChannels(sequence)
  maxComponents = maximum([length(channel.components) for channel in channels])
  result = zeros(Int64, (dfNumChannels(sequence), maxComponents))
  for (channelIdx, channel) in enumerate(channels)
    for (componentIdx, component) in enumerate(channel.components)
      result[channelIdx, componentIdx] = component.divider
    end
  end
  return result
end

dfNumChannels(sequence::Sequence) = length(periodicElectricalTxChannels(sequence)) # TODO: How do we integrate the mechanical channels?

function dfPhase(sequence::Sequence) # TODO: How do we integrate the mechanical channels and non-periodic channels and sweeps?
  channels = periodicElectricalTxChannels(sequence)
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

rxBandwidth(sequence::Sequence) = sequence.acquisiton.bandwidth
rxNumChannels(sequence::Sequence) = length(rxChannels(sequence))
rxNumSamplingPoints(sequence::Sequence) = round(Int64, upreferred(rxBandwidth(sequence)*2*dfCycle(sequence)))
rxNumSamplesPerPeriod(sequence::Sequence) = rxNumSamplingPoints(sequence)
rxChannels(sequence::Sequence) = sequence.acquisiton.channels

needsControl(sequence::Sequence) = any([field.control for field in sequence.fields])
needsDecoupling(sequence::Sequence) = any([field.decouple for field in sequence.fields])
needsControlOrDecoupling(sequence::Sequence) = needsControl(sequence) || needsDecoupling(sequence)