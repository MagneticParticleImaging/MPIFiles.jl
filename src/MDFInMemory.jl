using Unitful

export MDFv2InMemoryPart, MDFv2Variables, MDFv2InMemoryPart, MDFv2Root, MDFv2Study,
       MDFv2Experiment, MDFv2Tracer, MDFv2Scanner, MDFv2Drivefield, MDFv2Receiver,
       MDFv2Acquisition, MDFv2Measurement, MDFv2Calibration, MDFv2Reconstruction, MDFv2InMemory

export defaultMDFv2Root, defaultMDFv2Study, defaultMDFv2Experiment,
       defaultMDFv2Tracer, defaultMDFv2Scanner, defaultMDFv2Drivefield,
       defaultMDFv2Receiver, defaultMDFv2Acquisition, defaultMDFv2Measurement,
       defaultMDFv2Calibration, defaultMDFv2Reconstruction, defaultMDFv2InMemory

export checkConsistency, inMemoryMDFToDict, inMemoryMDFFromDict

abstract type MDFv2InMemoryPart end

"""
    $(TYPEDEF)

Set of variables of the MDF.
"""
mutable struct MDFv2Variables
  "tracer materials/injections for multi-color MPI"
  A::Union{Int64, Nothing}
  "acquired frames (N = O + E), same as a spatial position for calibration"
  N::Union{Int64, Nothing}
  "acquired background frames (E = N − O)"
  E::Union{Int64, Nothing}
  "acquired foreground frames (O = N − E)"
  O::Union{Int64, Nothing}
  "coefficients stored after sparsity transformation (B \\le O)"
  B::Union{Int64, Nothing}
  "periods within one frame"
  J::Union{Int64, Nothing}
  "partitions of each patch position"
  Y::Union{Int64, Nothing}
  "receive channels"
  C::Union{Int64, Nothing}
  "drive-field channels"
  D::Union{Int64, Nothing}
  "frequencies describing the drive-field waveform"
  F::Union{Int64, Nothing}
  "points sampled at receiver during one drive-field period"
  V::Union{Int64, Nothing}
  "sampling points containing processed data (W = V if no frequency selection or bandwidth reduction has been applied)"
  W::Union{Int64, Nothing}
  "frequencies describing the processed data (K = V/2 + 1 if no frequency selection or bandwidth reduction has been applied)"
  K::Union{Int64, Nothing}
  "frames in the reconstructed MPI data set"
  Q::Union{Int64, Nothing}
  "voxels in the reconstructed MPI data set"
  P::Union{Int64, Nothing}
  "channels in the reconstructed MPI data set"
  S::Union{Int64, Nothing}

  function MDFv2Variables()
    return new(nothing, nothing, nothing, nothing, nothing, nothing, nothing,
               nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
  end
end

"""
    $(TYPEDEF)
    
Root group of an in-memory MDF
"""
mutable struct MDFv2Root <: MDFv2InMemoryPart
  "UTC creation time of MDF data set"
  time::Union{DateTime, Missing}
  "Universally Unique Identifier (RFC 4122) of MDF file"
  uuid::Union{UUID, Missing}
  "Version of the file format"
  version::Union{VersionNumber, Missing}

  function MDFv2Root(;
    time = missing,
    uuid = missing,
    version = missing)
    
    return new(
      time,
      uuid,
      version
    )
  end
end

defaultMDFv2Root() = MDFv2Root(time=Dates.now(), uuid=UUIDs.uuid4(), version=VersionNumber("2.1.0"))

"""
    $(TYPEDEF)

Study group of an in-memory MDF
"""
mutable struct MDFv2Study <: MDFv2InMemoryPart
  "Short description of the study"
  description::Union{String, Missing}
  "Name of the study"
  name::Union{String, Missing}
  "Number of the study"
  number::Union{Int64, Missing}
  "UTC creation time of study; optional"
  time::Union{DateTime, Nothing}
  "Universally Unique Identifier (RFC 4122) of study"
  uuid::Union{UUID, Missing}

  function MDFv2Study(;
    description = missing,
    name = missing,
    number = missing,
    time = nothing,
    uuid = missing)
    
    return new(
      description,
      name,
      number,
      time,
      uuid
    )
  end
end

defaultMDFv2Study() = MDFv2Study(time=Dates.now(), uuid=UUIDs.uuid4())

"""
    $(TYPEDEF)

Experiment group of an in-memory MDF
"""
mutable struct MDFv2Experiment <: MDFv2InMemoryPart
  "Short description of the experiment"
  description::Union{String, Missing}
  "Flag indicating if the data in this file is simulated rather than measured"
  isSimulation::Union{Bool, Missing}
  "Experiment name"
  name::Union{String, Missing}
  "Experiment number within study"
  number::Union{Int64, Missing}
  "Name of the subject that was imaged"
  subject::Union{String, Missing}
  "Universally Unique Identifier (RFC 4122) of experiment"
  uuid::Union{UUID, Missing}

  function MDFv2Experiment(;
    description = missing,
    isSimulation = missing,
    name = missing,
    number = missing,
    subject = missing,
    uuid = missing)
    
    return new(
      description,
      isSimulation,
      name,
      number,
      subject,
      uuid
    )
  end
end

defaultMDFv2Experiment() = MDFv2Experiment(uuid=UUIDs.uuid4())

"""
    $(TYPEDEF)

Tracer group of an in-memory MDF; optional
"""
mutable struct MDFv2Tracer <: MDFv2InMemoryPart
  "Batch of tracer"
  batch::Union{Vector{String}, Missing}
  "A mol(solute)/L no Molar concentration of solute per litre"
  concentration::Union{Vector{Float64}, Missing}
  "UTC time at which tracer injection started; optional"
  injectionTime::Union{Vector{DateTime}, Nothing}
  "Name of tracer used in experiment"
  name::Union{Vector{String}, Missing}
  "Solute, e.g. Fe"
  solute::Union{Vector{String}, Missing}
  "Name of tracer supplier"
  vendor::Union{Vector{String}, Missing}
  "Total volume of applied tracer"
  volume::Union{Vector{Float64}, Missing}

  function MDFv2Tracer(;
    batch = missing,
    concentration = missing,
    injectionTime = nothing,
    name = missing,
    solute = missing,
    vendor = missing,
    volume = missing)
    
    return new(
      batch,
      concentration,
      injectionTime,
      name,
      solute,
      vendor,
      volume
    )
  end
end

defaultMDFv2Tracer() = MDFv2Tracer()

"""
    $(TYPEDEF)

Scanner group of an in-memory MDF
"""
mutable struct MDFv2Scanner <: MDFv2InMemoryPart
  "Diameter of the bore; optional"
  boreSize::Union{Float64, Nothing}
  "Facility where the MPI scanner is installed"
  facility::Union{String, Missing}
  "Scanner manufacturer"
  manufacturer::Union{String, Missing}
  "Scanner name"
  name::Union{String, Missing}
  "User who operates the MPI scanner"
  operator::Union{String, Missing}
  "Scanner topology (e.g. FFP, FFL, MPS)"
  topology::Union{String, Missing}

  function MDFv2Scanner(;
    boreSize = nothing,
    facility = missing,
    manufacturer = missing,
    name = missing,
    operator = missing,
    topology = missing)
    
    return new(
      boreSize,
      facility,
      manufacturer,
      name,
      operator,
      topology
    )
  end
end

defaultMDFv2Scanner() = MDFv2Scanner()

"""
    $(TYPEDEF)

Drivefield subgroup of acquisition group of an in-memory MDF
"""
mutable struct MDFv2Drivefield <: MDFv2InMemoryPart
  "Base frequency to derive drive field frequencies"
  baseFrequency::Union{Float64, Missing}
  "Trajectory cycle is determined by lcm(divider)/baseFrequency. It will not change
  when averaging was applied. The duration for measuring the V data points (i.e. the
  drive-field period) is given by the product of period and numAverages"
  cycle::Union{Float64, Missing}
  "Divider of the baseFrequency to determine the drive field frequencies"
  divider::Union{Array{Int64, 2}, Missing}
  "Number of drive field channels, denoted by D"
  numChannels::Union{Int64, Missing}
  "Applied drive field phase"
  phase::Union{Array{Float64, 3}, Missing}
  "Applied drive field strength"
  strength::Union{Array{Float64, 3}, Missing}
  "Waveform type: sine, triangle or custom"
  waveform::Union{Array{String, 2}, Missing}

  function MDFv2Drivefield(;
    baseFrequency = missing,
    cycle = missing,
    divider = missing,
    numChannels = missing,
    phase = missing,
    strength = missing,
    waveform = missing)
    
    return new(
      baseFrequency,
      cycle,
      divider,
      numChannels,
      phase,
      strength,
      waveform
    )
  end
end

defaultMDFv2Drivefield() = MDFv2Drivefield()

"""
    $(TYPEDEF)

Receiver subgroup of acquisition group of an in-memory MDF
"""
mutable struct MDFv2Receiver <: MDFv2InMemoryPart
  "Bandwidth of the receiver unit"
  bandwidth::Union{Float64, Missing}
  "Dimension less scaling factor and offset (a_c, b_c) to convert raw data into a
  physical quantity with corresponding unit of measurement unit; optional"
  dataConversionFactor::Union{Array{Float64, 2}, Nothing}
  "Induction factor mapping the projection of the magnetic moment to the voltage in the receive coil; optional"
  inductionFactor::Union{Vector{Float64}, Nothing}
  "Number of receive channels C"
  numChannels::Union{Int64, Missing}
  "Number of sampling points during one period, denoted by V"
  numSamplingPoints::Union{Int64, Missing}
  "Transfer function of the receive channels in Fourier domain. unit is the field
  from the /measurement group; optional"
  transferFunction::Union{Array{ComplexF64, 2}, Nothing}
  "SI unit of the measured quantity, usually Voltage V"
  unit::Union{String, Missing}

  function MDFv2Receiver(;
    bandwidth = missing,
    dataConversionFactor = nothing,
    inductionFactor = nothing,
    numChannels = missing,
    numSamplingPoints = missing,
    transferFunction = nothing,
    unit = missing)
    
    return new(
      bandwidth,
      dataConversionFactor,
      inductionFactor,
      numChannels,
      numSamplingPoints,
      transferFunction,
      unit
    )
  end
end

defaultMDFv2Receiver() = MDFv2Receiver(unit = "V")

"""
    $(TYPEDEF)

Acquisition group of an in-memory MDF
"""
mutable struct MDFv2Acquisition <: MDFv2InMemoryPart
  "Gradient strength of the selection field in x, y, and z directions; optional"
  gradient::Union{Array{Float64, 4}, Nothing}
  "Number of block averages per drive-field period."
  numAverages::Union{Int64, Missing}
  "Number of available measurement frames"
  numFrames::Union{Int64, Missing}
  "Number of drive-field periods within a frame denoted by J"
  numPeriodsPerFrame::Union{Int64, Missing}
  "Offset field applied; optional"
  offsetField::Union{Array{Float64, 3}, Nothing}
  "UTC start time of MPI measurement"
  startTime::Union{DateTime, Missing}

  drivefield::Union{MDFv2Drivefield, Missing}
  receiver::Union{MDFv2Receiver, Missing}

  function MDFv2Acquisition(;
    gradient = nothing,
    numAverages = missing,
    numFrames = missing,
    numPeriodsPerFrame = missing,
    offsetField = nothing,
    startTime = missing,
    drivefield = missing,
    receiver = missing)
    
    return new(
      gradient,
      numAverages,
      numFrames,
      numPeriodsPerFrame,
      offsetField,
      startTime,
      drivefield,
      receiver
    )
  end
end

defaultMDFv2Acquisition() = MDFv2Acquisition(startTime=Dates.now(), drivefield=defaultMDFv2Drivefield(), receiver=defaultMDFv2Receiver())

"""
    $(TYPEDEF)

Measurement group of an in-memory MDF
"""
mutable struct MDFv2Measurement <: MDFv2InMemoryPart
  "Measured data at a specific processing stage"
  data::Union{AbstractArray{<:Number, 4}, Missing}
  "Indices of original frame order; optional if !isFramePermutation"
  framePermutation::Union{Vector{Int64}, Nothing}
  "Indices of selected frequency components; optional if !isFrequencySelection"
  frequencySelection::Union{Vector{Int64}, Nothing}
  "Flag, if the background has been subtracted"
  isBackgroundCorrected::Union{Bool, Missing}
  "Mask indicating for each of the N frames if it is a background measurement (true) or not"
  isBackgroundFrame::Union{Vector{Bool}, Missing}
  "Flag, if the frame dimension N has been moved to the last dimension"
  isFastFrameAxis::Union{Bool, Missing}
  "Flag, if the data is stored in frequency space"
  isFourierTransformed::Union{Bool, Missing}
  "Flag, if the order of frames has been changed, see framePermutation"
  isFramePermutation::Union{Bool, Missing}
  "Flag, if only a subset of frequencies has been selected and stored, see frequencySelection"
  isFrequencySelection::Union{Bool, Missing}
  "Flag, if the foreground frames are compressed along the frame dimension"
  isSparsityTransformed::Union{Bool, Missing}
  "Flag, if spectral leakage correction has been applied"
  isSpectralLeakageCorrected::Union{Bool, Missing}
  "Flag, if the data has been corrected by the transferFunction"
  isTransferFunctionCorrected::Union{Bool, Missing}
  "Name of the applied sparsity transformation; optional if !isSparsityTransformed"
  sparsityTransformation::Union{String, Nothing}
  "Subsampling indices \\beta{j,c,k,b}; optional if !isSparsityTransformed"
  subsamplingIndices::Union{Array{Integer, 4}, Nothing}

  function MDFv2Measurement(;
    data = missing,
    framePermutation = nothing,
    frequencySelection = nothing,
    isBackgroundCorrected = missing,
    isBackgroundFrame = missing,
    isFastFrameAxis = missing,
    isFourierTransformed = missing,
    isFramePermutation = missing,
    isFrequencySelection = missing,
    isSparsityTransformed = missing,
    isSpectralLeakageCorrected = missing,
    isTransferFunctionCorrected = missing,
    sparsityTransformation = nothing,
    subsamplingIndices = nothing)

    return new(
      data,
      framePermutation,
      frequencySelection,
      isBackgroundCorrected,
      isBackgroundFrame,
      isFastFrameAxis,
      isFourierTransformed,
      isFramePermutation,
      isFrequencySelection,
      isSparsityTransformed,
      isSpectralLeakageCorrected,
      isTransferFunctionCorrected,
      sparsityTransformation,
      subsamplingIndices
    )
  end
end

defaultMDFv2Measurement() = MDFv2Measurement()

"""
    $(TYPEDEF)

Calibration group of an in-memory MDF
"""
mutable struct MDFv2Calibration <: MDFv2InMemoryPart
  "Size of the delta sample used for calibration scan; optional"
  deltaSampleSize::Union{Vector{Float64}, Nothing}
  "Field of view of the system matrix; optional"
  fieldOfView::Union{Vector{Float64}, Nothing}
  "Center of the system matrix (relative to origin/center); optional"
  fieldOfViewCenter::Union{Vector{Float64}, Nothing}
  "Method used to obtain calibration data. Can for instance be robot, hybrid or simulation"
  method::Union{String, Missing}
  "Applied offset field strength to emulate a spatial position (x, y, z); optional"
  offsetFields::Union{Array{Float64, 2}, Nothing}
  "Ordering of the dimensions, default is xyz; optional"
  order::Union{String, Nothing}
  "Position of each of the grid points, stored as (x, y, z) triples; optional"
  positions::Union{Array{Float64, 2}, Nothing}
  "Number of voxels in each dimension, inner product is O; optional"
  size::Union{Vector{Int64}, Nothing}
  "Signal-to-noise estimate for recorded frequency components; optional"
  snr::Union{Array{Float64, 3}, Nothing}
  "Flag, if the grid is meandering; optional"
  isMeanderingGrid::Union{Bool, Nothing}

  function MDFv2Calibration(;
    deltaSampleSize = nothing,
    fieldOfView = nothing,
    fieldOfViewCenter = nothing,
    method = missing,
    offsetFields = nothing,
    order = nothing,
    positions = nothing,
    size = nothing,
    snr = nothing,
    isMeanderingGrid = nothing)

    return new(
      deltaSampleSize,
      fieldOfView,
      fieldOfViewCenter,
      method,
      offsetFields,
      order,
      positions,
      size,
      snr,
      isMeanderingGrid
    )
  end
end

defaultMDFv2Calibration() = MDFv2Calibration(order="xyz")

"""
    $(TYPEDEF)

Reconstruction group of an in-memory MDF
"""
mutable struct MDFv2Reconstruction <: MDFv2InMemoryPart
  "Reconstructed data"
  data::Union{Array{Float32, 3}, Missing}
  "Field of view of reconstructed data; optional"
  fieldOfView::Union{Vector{Float64}, Nothing}
  "Center of the reconstructed data (relative to scanner origin/center); optional"
  fieldOfViewCenter::Union{Vector{Float64}, Nothing}
  "Mask indicating for each of the P voxels if it is part of the overscan region (true) or not; optional"
  isOverscanRegion::Union{Vector{Bool}, Nothing}
  "Ordering of the dimensions, default is xyz; optional"
  order::Union{String, Nothing}
  "Position of each of the grid points, stored as (x, y, z) tripels; optional"
  positions::Union{Array{Float64, 2}, Nothing}
  "Number of voxels in each dimension, inner product is P; optional"
  size::Union{Vector{Int64}, Nothing}

  function MDFv2Reconstruction(;
    data = missing,
    fieldOfView = nothing,
    fieldOfViewCenter = nothing,
    isOverscanRegion = nothing,
    order = nothing,
    positions = nothing,
    size = nothing)

    return new(
      data,
      fieldOfView,
      fieldOfViewCenter,
      isOverscanRegion,
      order,
      positions,
      size
    )
  end
end

defaultMDFv2Reconstruction() = MDFv2Reconstruction(order="xyz")

"""
    $(TYPEDEF)

In-memory description of an MDF file.
"""
mutable struct MDFv2InMemory <: MPIFile # TODO: Not sure, if MPIFile is a good fit
  root::Union{MDFv2Root, Missing}
  study::Union{MDFv2Study, Missing}
  experiment::Union{MDFv2Experiment, Missing}
  tracer::Union{MDFv2Tracer, Nothing}
  scanner::Union{MDFv2Scanner, Missing}
  acquisition::Union{MDFv2Acquisition, Missing}
  measurement::Union{MDFv2Measurement, Nothing}
  calibration::Union{MDFv2Calibration, Nothing}
  reconstruction::Union{MDFv2Reconstruction, Nothing}
  custom::Dict{String, Any}
  variables::MDFv2Variables

  function MDFv2InMemory(;
    root = missing,
    study = missing,
    experiment = missing,
    tracer = nothing,
    scanner = missing,
    acquisition = missing,
    measurement = nothing,
    calibration = nothing,
    reconstruction = nothing,
    custom = Dict{String, Any}(),
    variables = MDFv2Variables())

    return new(
      root,
      study,
      experiment,
      tracer,
      scanner,
      acquisition,
      measurement,
      calibration,
      reconstruction,
      custom,
      variables
    )
  end
end

MDFv2InMemory(mdfFile::MDFFileV2) = inMemoryMDFFromMDFFileV2(mdfFile)
MDFv2InMemory(dict::Dict) = inMemoryMDFFromDict(dict)

function defaultMDFv2InMemory()
  return MDFv2InMemory(
    root = defaultMDFv2Root(),
    study = defaultMDFv2Study(),
    experiment = defaultMDFv2Experiment(),
    scanner = defaultMDFv2Scanner(),
    acquisition = defaultMDFv2Acquisition()
  )
end

function check(part::MDFv2Root, variables::MDFv2Variables)
  @assert part.version > VersionNumber("2") "The version should be at least version 2."
end

function check(part::MDFv2Study, variables::MDFv2Variables)
  # nothing to check yet
end

function check(part::MDFv2Experiment, variables::MDFv2Variables)
  # nothing to check yet
end

function check(part::MDFv2Tracer, variables::MDFv2Variables)
  # Init comparison
  if isnothing(variables.A)
    variables.A = length(part.vendor)
  end

  # Check if all sizes match
  for fieldname in fieldnames(MDFv2Tracer)
    field = getproperty(part, fieldname)
    if !isnothing(field)
      @assert length(field) == variables.A "Inconsistent dimensions in `$fieldname` in `tracer`."
    end
  end
end

function check(part::MDFv2Scanner, variables::MDFv2Variables)
  # nothing to check yet
end

function check(part::MDFv2Acquisition, variables::MDFv2Variables)
  # Pick variables first
  if isnothing(variables.N)
    variables.N = part.numFrames
  end
  if isnothing(variables.J)
    variables.J = part.numPeriodsPerFrame
  end

  # Check dimensions of `gradient` field
  if !isnothing(part.gradient)
    if isnothing(variables.J)
      variables.J = size(part.gradient, 4)
    else
      @assert variables.J == size(part.gradient, 4) "Inconsistent dimension J in `gradient` in `acquisition`."
    end

    if isnothing(variables.Y)
      variables.Y = size(part.gradient, 3)
    else
      @assert variables.Y == size(part.gradient, 3) "Inconsistent dimension Y in `gradient` in `acquisition`."
    end

    @assert size(part.gradient, 1) == 3 "Inconsistent first dimension in `gradient` in `acquisition`."
    @assert size(part.gradient, 2) == 3 "Inconsistent second dimension in `gradient` in `acquisition`."
  end

  # Check dimensions of `offsetField` field
  if !isnothing(part.offsetField)
    if isnothing(variables.J)
      variables.J = size(part.offsetField, 3)
    else
      @assert variables.J == size(part.offsetField, 3) "Inconsistent dimension J in `offsetField` in `acquisition`."
    end

    if isnothing(variables.Y)
      variables.Y = size(part.offsetField, 2)
    else
      @assert variables.Y == size(part.offsetField, 2) "Inconsistent dimension Y in `offsetField` in `acquisition`."
    end

    @assert size(part.offsetField, 1) == 3 "Inconsistent first dimension in `offsetField` in `acquisition`."
  end
end

function check(part::MDFv2Drivefield, variables::MDFv2Variables)
  # Pick variables first
  if isnothing(variables.D)
    variables.D = part.numChannels
  end
  if isnothing(variables.F)
    variables.F = size(part.divider, 1)
  end
  if isnothing(variables.J)
    variables.J = size(part.phase, 3)
  end

  # Then check dimensions for multidimensional fields
  @assert variables.F == size(part.divider, 1) "Inconsistent dimension F in `divider` in `drivefield`."
  @assert variables.D == size(part.divider, 2) "Inconsistent dimension D in `divider` in `drivefield`."

  @assert variables.F == size(part.phase, 1) "Inconsistent dimension F in `phase` in `drivefield`."
  @assert variables.D == size(part.phase, 2) "Inconsistent dimension D in `phase` in `drivefield`."
  @assert variables.J == size(part.phase, 3) "Inconsistent dimension J in `phase` in `drivefield`."

  @assert variables.F == size(part.strength, 1) "Inconsistent dimension F in `strength` in `drivefield`."
  @assert variables.D == size(part.strength, 2) "Inconsistent dimension D in `strength` in `drivefield`."
  @assert variables.J == size(part.strength, 3) "Inconsistent dimension J in `strength` in `drivefield`."

  @assert variables.F == size(part.waveform, 1) "Inconsistent dimension F in `waveform` in `drivefield`."
  @assert variables.D == size(part.waveform, 2) "Inconsistent dimension D in `waveform` in `drivefield`."
end

function check(part::MDFv2Receiver, variables::MDFv2Variables)
  # C is defined by numChannels
  if isnothing(variables.C)
    variables.C = part.numChannels
  end
  # V is defined by numSamplingPoints
  if isnothing(variables.V)
    variables.V = part.numSamplingPoints
  end

  if !isnothing(part.dataConversionFactor)
    @assert variables.C == size(part.dataConversionFactor, 2) "Inconsistent dimension C in `dataConversionFactor` in `receiver`."
    @assert size(part.dataConversionFactor, 1) == 2 "Inconsistent first dimension in `dataConversionFactor` in `receiver`."
  end

  if !isnothing(part.inductionFactor)
    @assert variables.C == length(part.inductionFactor) "Inconsistent dimension C in `inductionFactor` in `receiver`."
  end

  if !isnothing(part.transferFunction)
    @assert variables.C == size(part.transferFunction, 2) "Inconsistent dimension C in `transferFunction` in `receiver`."

    if isnothing(variables.K)
      variables.K = size(part.transferFunction, 1)
    else
      @assert variables.K == size(part.transferFunction, 1) "Inconsistent dimension K in `transferFunction` in `receiver`."
    end
  end
end

function check(part::MDFv2Measurement, variables::MDFv2Variables)
  # TODO: check data consistency; I don't know how to do that since we can't know all variables
  # N × J × C × K or
  # J × C × K × N or
  # N × J × C × W or
  # J × C × W × N or
  # J×C×K×(B+E)

  # N is defined by the length of isBackgroundFrame here; should have been retrieved earlier
  if isnothing(variables.N)
    variables.N = length(part.isBackgroundFrame)
  end
  # E is defined by the number of `true` values in isBackgroundFrame here; should have been retrieved earlier
  if isnothing(variables.E)
    variables.E = length([x for x in part.isBackgroundFrame if x == true])
  end
  # E is defined by the number of `false` values in isBackgroundFrame here; should have been retrieved earlier
  if isnothing(variables.O)
    variables.O = length([x for x in part.isBackgroundFrame if x == false])
  end
  if isnothing(variables.W)
    if isnothing(variables.V)
      @warn "Can't determine variable W since variable V is not defined. Should happen in receiver."
    else
      variables.W = variables.V
    end
  end

  if isnothing(variables.J)
    if !part.isSparsityTransformed
      @debug "Can't determine variable J for measurement from `subsamplingIndices` since the dataset is not sparsity transformed."
    else
      variables.J = size(part.subsamplingIndices, 1)
    end
  end
  if isnothing(variables.C)
    if !part.isSparsityTransformed
      @debug "Can't determine variable C for measurement from `subsamplingIndices` since the dataset is not sparsity transformed. Should happen in receiver."
    else
      variables.C = size(part.subsamplingIndices, 2)
    end
  end
  if isnothing(variables.K)
    if !part.isSparsityTransformed
      @debug "Can't determine variable K for measurement from `subsamplingIndices` since the dataset is not sparsity transformed. Should happen in receiver."
    else
      variables.K = size(part.subsamplingIndices, 3)
    end
  end
  if isnothing(variables.B)
    if !part.isSparsityTransformed
      @debug "Can't determine variable B for measurement from `subsamplingIndices` since the dataset is not sparsity transformed. Should happen in receiver."
    else
      variables.B = size(part.subsamplingIndices, 4)
    end
  end


  if part.isFramePermutation
    @assert variables.N == length(part.framePermutation) "Inconsistent dimension N in `framePermutation` in `measurement`."
  end

  if part.isFrequencySelection
    if isnothing(variables.K)
      variables.K = length(part.frequencySelection)
    else
      @assert variables.K == length(part.frequencySelection) "Inconsistent dimension K in `frequencySelection` in `measurement`."
    end
  end

  @assert variables.N == length(part.isBackgroundFrame) "Inconsistent dimension N in `isBackgroundFrame` in `measurement`."

  if part.isSparsityTransformed
    @assert !isnothing(part.sparsityTransformation) "Field `sparsityTransformation` must be set when `isSparsityTransformed` is set in in `measurement`."
  end

  if part.isSparsityTransformed
    # J, C, K and B should be defined by now
    @assert variables.B == size(part.subsamplingIndices, 1) "Inconsistent dimension B in `subsamplingIndices` in `measurement`."
    @assert variables.K == size(part.subsamplingIndices, 2) "Inconsistent dimension K in `subsamplingIndices` in `measurement`."
    @assert variables.C == size(part.subsamplingIndices, 3) "Inconsistent dimension C in `subsamplingIndices` in `measurement`."
    @assert variables.J == size(part.subsamplingIndices, 4) "Inconsistent dimension J in `subsamplingIndices` in `measurement`."
  end
end

function check(part::MDFv2Calibration, variables::MDFv2Variables)
  if !isnothing(part.deltaSampleSize)
    @assert length(part.deltaSampleSize) == 3 "Inconsistent length in `deltaSampleSize` in `calibration`."
  end

  if !isnothing(part.fieldOfView)
    @assert length(part.fieldOfView) == 3 "Inconsistent length in `fieldOfView` in `calibration`."
  end

  if !isnothing(part.fieldOfViewCenter)
    @assert length(part.fieldOfViewCenter) == 3 "Inconsistent length in `fieldOfViewCenter` in `calibration`."
  end

  if !isnothing(part.offsetFields)
    if isnothing(variables.O)
      variables.O = size(part.offsetFields, 2)
    else
      @assert variables.O == size(part.offsetFields, 2) "Inconsistent dimension O in `offsetFields` in `calibration`."
    end
    @assert size(part.offsetFields, 1) == 3 "Inconsistent first dimension in `offsetFields` in `calibration`."
  end

  if !isnothing(part.order)
    @assert part.order in ["xyz", "xzy", "yxz", "yzx", "zyx", "zxy"] "Wrong `order` of `$(part.order)` in `calibration`."
  end

  if !isnothing(part.positions)
    # O must be defined by now
    @assert variables.O == size(part.positions, 2) "Inconsistent dimension O in `positions` in `calibration`."
    @assert size(part.positions, 1) == 3 "Inconsistent first dimension in `positions` in `calibration`."
  end

  if !isnothing(part.size)
    @assert length(part.size) == 3 "Inconsistent length in `size` in `calibration`."
    # O must be defined by now
    @assert variables.O == prod(part.size) "The product of `size` with `$(part.size)` must equal O."
  end

  if !isnothing(part.snr)
    if isnothing(variables.J)
      variables.J = size(part.snr, 3)
    else
      @assert variables.J == size(part.snr, 3) "Inconsistent dimension J in `snr` in `calibration`."
    end

    if isnothing(variables.C)
      variables.C = size(part.snr, 2)
    else
      @assert variables.C == size(part.snr, 2) "Inconsistent dimension C in `snr` in `calibration`."
    end

    if isnothing(variables.K)
      variables.K = size(part.snr, 1)
    else
      @assert variables.K == size(part.snr, 1) "Inconsistent dimension K in `snr` in `calibration`."
    end
  end
end

function check(part::MDFv2Reconstruction, variables::MDFv2Variables)
  # Pick variables first
  if isnothing(variables.Q)
    variables.Q = size(part.data, 3)
  end
  if isnothing(variables.P)
    variables.P = size(part.data, 2)
  end
  if isnothing(variables.S)
    variables.S = size(part.data, 1)
  end
  
  if !isnothing(part.fieldOfView)
    @assert length(part.fieldOfView) == 3 "Inconsistent length in `fieldOfView` in `reconstruction`."
  end

  if !isnothing(part.fieldOfViewCenter)
    @assert length(part.fieldOfViewCenter) == 3 "Inconsistent length in `fieldOfViewCenter` in `reconstruction`."
  end

  if !isnothing(part.isOverscanRegion)
    @assert variables.P == length(part.isOverscanRegion) "Inconsistent length in `isOverscanRegion` in `reconstruction`."
  end

  if !isnothing(part.order)
    @assert part.order in ["xyz", "xzy", "yxz", "yzx", "zyx", "zxy"] "Wrong `order` of `$(part.order)` in `reconstruction`."
  end

  if !isnothing(part.positions)
    @assert variables.P == size(part.positions, 2) "Inconsistent dimension P in `positions` in `reconstruction`."
    @assert size(part.positions, 1) == 3 "Inconsistent first dimension in `positions` in `reconstruction`."
  end

  if !isnothing(part.size)
    @assert length(part.size) == 3 "Inconsistent length in `size` in `reconstruction`."
    # P must be defined by now
    @assert variables.P == prod(part.size) "The product of `size` with `$(part.size)` must equal P."
  end
end

"Check for missing fields in MDF parts"
function checkMissing(part::T) where T <: MDFv2InMemoryPart
  for (fieldname, fieldtype) in zip(fieldnames(T), fieldtypes(T))
    if !(fieldtype == MDFv2Variables)
      field = getproperty(part, fieldname)
      @assert !ismissing(field) "The field `$fieldname` is missing in the given in-memory MDF."
    end
  end
end

"Check, whether all non-optional fields have been
set and if the dimensions of the fields match"
function checkConsistency(mdf::MDFv2InMemory)
  for (fieldname, fieldtype) in zip(fieldnames(MDFv2InMemory), fieldtypes(MDFv2InMemory))
    if !(fieldtype == MDFv2Variables || fieldname == :custom)
      # At the moment, this should be a Union
      fieldtype = (fieldtype.b <: MDFv2InMemoryPart) ? fieldtype.b : fieldtype.a

      @debug "Checking consistency of `$fieldname`."

      field = getproperty(mdf, fieldname)
      @assert !ismissing(field) "The field `$fieldname` is missing in the given in-memory MDF."
       
      if !isnothing(field) # Only optional fields can be `nothing`
        checkMissing(field)
        check(field, mdf.variables)
      end

      # Check subgroups of acquisition
      if fieldtype == MDFv2Acquisition
        checkMissing(field.drivefield)
        check(field.drivefield, mdf.variables)

        checkMissing(field.receiver)
        check(field.receiver, mdf.variables)
      end
    end
  end

  return true # If no assertions failed, return true for checking in tests
end

### Create getters and setters

customSymbols = Dict{Symbol, String}(
  :dfCustomWaveform => "/acquisition/drivefield/customWaveform",
  :measTemperatures => "/measurement/_monitoring/temperature/observed",
  :measObservedDriveField => "/measurement/_monitoring/driveField/observed",
  :measAppliedDriveField => "/measurement/_monitoring/driveField/applied",
  :rxTransferFunctionFileName => "/acquisition/receiver/transferFunctionFileName",
  :recoParameters => "/reconstruction/_parameters",
  :auxiliaryData => "/custom/auxiliaryData",
)

specificationSymbols = Dict{Symbol, String}()

aliases = Dict{Symbol, Symbol}(
  :measIsTransferFunctionCorrected => :measIsTFCorrected,
  :measIsBackgroundCorrected => :measIsBGCorrected,
  :measIsBackgroundFrame => :measIsBGFrame,
  :calibSnr => :calibSNR,
  :calibFieldOfView => :calibFov,
  :calibFieldOfViewCenter => :calibFovCenter,
  :recoFieldOfView => :recoFov,
  :recoFieldOfViewCenter => :recoFovCenter
)

prefixes = Dict{String, String}(
  "MDFv2Root" => "",
  "MDFv2Study" => "study",
  "MDFv2Experiment" => "experiment",
  "MDFv2Tracer" => "tracer",
  "MDFv2Scanner" => "scanner",
  "MDFv2Acquisition" => "acq",
  "MDFv2Drivefield" => "df",
  "MDFv2Receiver" => "rx",
  "MDFv2Measurement" => "meas",
  "MDFv2Calibration" => "calib",
  "MDFv2Reconstruction" => "reco"
)
for (fieldname, fieldtype) in zip(fieldnames(MDFv2InMemory), fieldtypes(MDFv2InMemory))
  fieldnameStr = string(fieldname)
  if !(fieldnameStr == "variables" || fieldnameStr == "custom")
    # At the moment, this should be a Union
    missingOrNothing = (fieldtype.b <: MDFv2InMemoryPart) ? fieldtype.a : fieldtype.b
    fieldtype = (fieldtype.b <: MDFv2InMemoryPart) ? fieldtype.b : fieldtype.a

    capitalizedFieldnameStr = uppercase(fieldnameStr[1])*fieldnameStr[2:end]

    # Create getter and setter for the whole group
    @eval begin
      export $fieldname
      @doc $"""
          $fieldnameStr(mdf)

      $capitalizedFieldnameStr group of an in-memory MDF.
      """
      function $(fieldname)(mdf::MDFv2InMemory)::Union{$fieldtype, $missingOrNothing}
        return mdf.$fieldname
      end

      @doc $"""
          $fieldnameStr(mdf, value)

      $capitalizedFieldnameStr group of an in-memory MDF.
      """
      function $(fieldname)(mdf::MDFv2InMemory, value::Union{$fieldtype, $missingOrNothing})
        mdf.$fieldname = value
      end
    end

    for (partFieldname, partFieldtype) in zip(fieldnames(fieldtype), fieldtypes(fieldtype))
      partFieldnameStr = string(partFieldname)

      fieldDocstring = fielddoc(fieldtype, partFieldname)

      # The acquisition group has subgroups, so we need to go deeper there
      if !(partFieldnameStr == "drivefield" || partFieldnameStr == "receiver")
        if fieldtype != MDFv2Root
          capitalizedPartFieldname = uppercase(partFieldnameStr[1])*partFieldnameStr[2:end]
        else
          capitalizedPartFieldname = partFieldnameStr
        end
        functionSymbol = Symbol(prefixes[replace(string(fieldtype), "MPIFiles." => "")]*capitalizedPartFieldname)

        # Save symbols for later use in conversion
        if fieldnameStr != "root"
          specificationSymbols[functionSymbol] = "/"*fieldnameStr*"/"*partFieldnameStr
        else
          specificationSymbols[functionSymbol] = "/"*partFieldnameStr
        end

        @eval begin
          @doc $"""
              $functionSymbol(mdfPart)

          $fieldDocstring
          """
          function $(functionSymbol)(mdfPart::$fieldtype)::Union{$partFieldtype, $missingOrNothing}
            return mdfPart.$partFieldname
          end

          @doc $"""
              $functionSymbol(mdfPart, value)

          $fieldDocstring
          """
          function $(functionSymbol)(mdfPart::$fieldtype, value::Union{$partFieldtype, $missingOrNothing})
            mdfPart.$partFieldname = value
          end

          @doc $"""
              $functionSymbol(mdf)

          $fieldDocstring
          """
          function $(functionSymbol)(mdf::MDFv2InMemory)::Union{$partFieldtype, $missingOrNothing}
            if !(isnothing($fieldname(mdf)) || ismissing($fieldname(mdf)))
              return $(functionSymbol)($fieldname(mdf))
            else
              return $fieldname(mdf)
            end
          end

          @doc $"""
              $functionSymbol(mdf, value)

          $fieldDocstring
          """
          function $(functionSymbol)(mdf::MDFv2InMemory, value::Union{$partFieldtype, $missingOrNothing})
            # Automatically create fields if they do not exist
            if isnothing($fieldname(mdf)) || ismissing($fieldname(mdf))
              @debug "Creating field $($fieldnameStr)"
              $fieldname(mdf, $fieldtype())
            end

            $(functionSymbol)($fieldname(mdf), value)
          end
        end

        # If needed, create aliases
        if haskey(aliases, functionSymbol)
          alias = aliases[functionSymbol]
          @eval begin
            @doc $"""
              $alias(mdf)

            $fieldDocstring
            """
            $(alias)(mdf::MDFv2InMemory)::Union{$partFieldtype, $missingOrNothing} = $(functionSymbol)(mdf)

            @doc $"""
              $alias(mdf, value)

            $fieldDocstring
            """
            $(alias)(mdf::MDFv2InMemory, value::$partFieldtype) = $(functionSymbol)(mdf, value)

            $(functionSymbol)(f::MDFFileV2)::Union{$partFieldtype, $missingOrNothing} = $(alias)(f) # Should this be here or in MDF.jl?
          end
        end
      else
        # At the moment, this should be a Union
        missingOrNothing = (partFieldtype.b <: MDFv2InMemoryPart) ? partFieldtype.a : partFieldtype.b
        partFieldtype = (partFieldtype.b <: MDFv2InMemoryPart) ? partFieldtype.b : partFieldtype.a

        # Create getter and for the whole group
        @eval begin
          @doc $"""
              $partFieldname(mdfPart)

          $fieldDocstring
          """
          function $(partFieldname)(mdfPart::$fieldtype)::Union{$partFieldtype, $missingOrNothing}
            return mdfPart.$partFieldname
          end

          @doc $"""
              $partFieldname(mdfPart, value)

          $fieldDocstring
          """
          function $(partFieldname)(mdfPart::$fieldtype, value::Union{$partFieldtype, $missingOrNothing})
            mdfPart.$partFieldname = value
          end

          @doc $"""
              $partFieldname(mdf)

          $fieldDocstring
          """
          function $(partFieldname)(mdf::MDFv2InMemory)::Union{$partFieldtype, $missingOrNothing}
            if !(isnothing($fieldname(mdf)) || ismissing($fieldname(mdf)))
              return $partFieldname($fieldname(mdf))
            else
              return $fieldname(mdf)
            end
          end

          @doc $"""
              $partFieldname(mdfPart, value)

          $fieldDocstring
          """
          function $(partFieldname)(mdf::MDFv2InMemory, value::Union{$partFieldtype, $missingOrNothing})
            # Automatically create fields if they do not exist
            if isnothing($fieldname(mdf)) || ismissing($fieldname(mdf))
              @debug "Creating field $($fieldnameStr)"
              $fieldname(mdf, $fieldtype())
            end

            $partFieldname($fieldname(mdf), value)
          end
        end

        for (subPartFieldname, subPartFieldtype) in zip(fieldnames(partFieldtype), fieldtypes(partFieldtype))
          subPartFieldnameStr = string(subPartFieldname)
          capitalizedSubPartFieldname = uppercase(subPartFieldnameStr[1])*subPartFieldnameStr[2:end]
          functionSymbol = Symbol(prefixes[replace(string(partFieldtype), "MPIFiles." => "")]*capitalizedSubPartFieldname)

          subPartMissingOrNothing = missingOrNothing # TODO: This should be type-specific

          # Save symbols for later use in conversion
          specificationSymbols[functionSymbol] = "/"*fieldnameStr*"/"*partFieldnameStr*"/"*subPartFieldnameStr

          subPartFieldDocstring = fielddoc(partFieldtype, subPartFieldname)

          # Create getter and setter for the respective field  within the naming scheme
          @eval begin
            @doc $"""
              $functionSymbol(mdfPart)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdfPart::$partFieldtype)::Union{$subPartFieldtype, $subPartMissingOrNothing}
              return mdfPart.$subPartFieldname
            end

            @doc $"""
              $functionSymbol(mdfPart, value)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdfPart::$partFieldtype, value::Union{$subPartFieldtype, $subPartMissingOrNothing})
              mdfPart.$subPartFieldname = value
            end

            @doc $"""
              $functionSymbol(mdfPart)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdfPart::$fieldtype)::Union{$partFieldtype, $subPartMissingOrNothing}
              if !(isnothing($partFieldname(mdfPart)) || ismissing($partFieldname(mdfPart)))
                return $subPartFieldname($partFieldname(mdfPart))
              else
                return $partFieldname(mdfPart)
              end
            end

            @doc $"""
              $functionSymbol(mdfPart, value)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdfPart::$fieldtype, value::Union{$partFieldtype, $subPartMissingOrNothing})
              # Automatically create fields if they do not exist
              if isnothing($partFieldname(mdfPart)) || ismissing($subPartFieldname(mdfPart))
                @debug "Creating field $($partFieldnameStr).$($subPartFieldnameStr)"
                $partFieldname(mdfPart, $partFieldtype())
              end
              
              $subPartFieldname(mdfPart, value)
            end

            @doc $"""
              $functionSymbol(mdf)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdf::MDFv2InMemory)::Union{$subPartFieldtype, $subPartMissingOrNothing}
              if !(isnothing($fieldname(mdf)) || ismissing($fieldname(mdf)))
                if !(isnothing($partFieldname($fieldname(mdf))) || ismissing($partFieldname($fieldname(mdf))))
                  return $functionSymbol($partFieldname($fieldname(mdf)))
                else
                  return $partFieldname($fieldname(mdf))
                end
              else
                return $fieldname(mdf)
              end
            end
  
            @doc $"""
              $functionSymbol(mdf, value)

            $subPartFieldDocstring
            """
            function $(functionSymbol)(mdf::MDFv2InMemory, value::Union{$subPartFieldtype, $subPartMissingOrNothing})
              # Automatically create fields if they do not exist
              if isnothing($fieldname(mdf)) || ismissing($fieldname(mdf))
                @debug "Creating field $($fieldnameStr)"
                $fieldname(mdf, $fieldtype())
              end

              if isnothing($partFieldname($fieldname(mdf))) || ismissing($partFieldname($fieldname(mdf)))
                @debug "Creating field $($fieldnameStr).$($partFieldnameStr)"
                $partFieldname($fieldname(mdf), $partFieldtype())
              end
              
              $functionSymbol($partFieldname($fieldname(mdf)), value)
            end
          end

          # If needed, create aliases
          if haskey(aliases, functionSymbol)
            alias = aliases[functionSymbol]
            @eval begin
              $(alias)(mdf::MDFv2InMemory)::Union{$partFieldtype, $missingOrNothing} = $(functionSymbol)(mdf)
              $(alias)(mdf::MDFv2InMemory, value::Union{$partFieldtype, $missingOrNothing}) = $(functionSymbol)(mdf, value)
              $(functionSymbol)(f::MDFFileV2)::Union{$partFieldtype, $missingOrNothing} = $(alias)(f) # Should this be here or in MDF.jl?
            end
          end
        end
      end
    end
  end
end

# And some utility functions
measIsCalibProcessed(mdf::MDFv2InMemory)::Union{Bool, Missing} = measIsFramePermutation(mdf) && 
                                                                 measIsFourierTransformed(mdf) &&
                                                                 measIsFastFrameAxis(mdf)

experimentHasReconstruction(mdf::MDFv2InMemory)::Bool = !isnothing(mdf.reconstruction)
experimentHasMeasurement(mdf::MDFv2InMemory)::Bool = !isnothing(mdf.measurement)
rxHasTransferFunction(mdf::MDFv2InMemory) = !isnothing(mdf.acquisition.receiver.transferFunction)

# Last part is a little workaround for the testcases, since the file somehow contains a `deltaSampleSize` field
experimentIsCalibration(mdf::MDFv2InMemory)::Bool = !isnothing(mdf.calibration) && !ismissing(mdf.calibration.method)

# Creation and conversion

"Create an in-memory MDF from a dict matching the respective function names."
function inMemoryMDFFromDict(dict::Dict{Symbol, Any})::MDFv2InMemory
  mdf = MDFv2InMemory()

  for (functionSymbol, value) in dict
    # Conversion uses aliases => replace these symbols with their alias
    if haskey(aliases, functionSymbol)
      functionSymbol = aliases[functionSymbol]
    end

    try
      f = getfield(MPIFiles, functionSymbol)
      f(mdf, value)
    catch e
      if e isa UndefVarError
        @warn "The function `$(string(functionSymbol))` corresponding to a key in the given dict could not be found."
      else
        rethrow()
      end
    end
  end

  return mdf
end
inMemoryMDFFromDict(dict::Dict{String, Any})::MDFv2InMemory = inMemoryMDFFromDict(Dict(map(x -> Symbol(x.first)=>x.second, collect(dict))))

"Create a dict from an in-memory MDF by matching the respective function names."
function inMemoryMDFToDict(mdf::MDFv2InMemory)::Dict{String, Any}
  resultDict = Dict{String, Any}()

  # Add standard-compliant and non-standard fields
  for functionSymbol in vcat(collect(keys(specificationSymbols)), collect(keys(customSymbols)))
    f = getfield(MPIFiles, functionSymbol)
    result = f(mdf)
    if !(isnothing(result) || ismissing(result))
      # Conversion uses aliases => replace these symbols with their alias
      if haskey(aliases, functionSymbol)
        functionSymbol = aliases[functionSymbol]
      end

      resultDict[string(functionSymbol)] = result
    end
  end

  # Add measurements data
  resultDict["measData"] = measDataRaw(mdf)

  return resultDict
end

"Alias to inMemoryMDFToDict()"
toDict(mdf::MDFv2InMemory)::Dict{String, Any} = inMemoryMDFToDict(mdf)

"Create an in-memory MDF from an MDFFile by calling the corresponding functions."
function inMemoryMDFFromMDFFileV2(mdfFile::MDFFileV2)::MDFv2InMemory
  mdf = MDFv2InMemory()

  # Add standard-compliant fields
  for functionSymbol in keys(specificationSymbols)
    f = getfield(MPIFiles, functionSymbol)

    result = nothing
    try
      result = f(mdfFile)
    catch e
      @warn "Exception while reading symbol $(functionSymbol). Please check closely. Exception was `$e`."
    end

    if !(isnothing(result) || ismissing(result))
      f(mdf, result) # Call the setter of an MDFv2InMemory with a getter from an MDFFileV2
    end
  end

  # Add non-standard fields
  for functionSymbol in keys(customSymbols)
    f = getfield(MPIFiles, functionSymbol)

    result = nothing
    try
      result = f(mdfFile)
    catch e
      @warn "Exception while reading symbol $(functionSymbol). Please check closely."
    end

    if !(isnothing(result) || ismissing(result))
      mdf.custom[string(functionSymbol)] = result
    end
  end

  # Add measurements data
  if !isnothing(mdfFile.mmap_measData)
    measDataRaw(mdf, mdfFile.mmap_measData)
  else
    @warn "The measurement data could not be read. Please check closely."
  end

  return mdf
end

"Create an MDFFile from an in-memory MDF; alias to `saveasMDF`."
function inMemoryMDFToMDFFileV2(filename::String, mdf::MDFv2InMemory)
  saveasMDF(filename, mdf)
end

function saveasMDF(filename::String, mdf::MDFv2InMemory; failOnInconsistent::Bool=false)
  # file has to be removed if exists. Otherwise h5create fails.
  isfile(filename) && rm(filename)
  h5open(filename, "w") do file
    saveasMDF(file, mdf)
  end
  return filename
end

"Create an MDFFile from an in-memory MDF."
function saveasMDF(file::HDF5.File, mdf::MDFv2InMemory; failOnInconsistent::Bool=false)
  try
    checkConsistency(mdf)
  catch e
    if e isa AssertionError && !failOnInconsistent
      @warn "There is an inconsistency in the given in-memory MDF. The message is: `$(e.msg)`."
    else
      rethrow()
    end
  end

  for (functionSymbol, key) in merge(customSymbols, specificationSymbols)
    f = getfield(MPIFiles, functionSymbol)
    result = f(mdf)
    if !(isnothing(result) || ismissing(result))
      # Convert datatypes to something compatible with HDF5
      if eltype(result) == Bool
        result = Int8.(result)
      end
      if result isa VersionNumber || result isa DateTime || result isa UUID
        result = string(result)
      end
      if eltype(result) == DateTime
        result = string.(result)
      end

      file[key] = result
    end
  end
end

# This is non-standard (add new non-standard functions to `customSymbols` in order to have the custom fields set!)
dfCustomWaveform(mdf::MDFv2InMemory)::Union{String, Nothing} = @keyoptional mdf.custom["dfCustomWaveform"] # TODO: Should this be a 2D array?
dfCustomWaveform(mdf::MDFv2InMemory, customWaveform::String) = mdf.custom["dfCustomWaveform"] = customWaveform

rxTransferFunctionFileName(mdf::MDFv2InMemory)::Union{String, Nothing} = @keyoptional mdf.custom["rxTransferFunctionFileName"]
rxTransferFunctionFileName(mdf::MDFv2InMemory, filename::String) = mdf.custom["rxTransferFunctionFileName"] = filename

recoParameters(mdf::MDFv2InMemory) = @keyoptional mdf.custom["recoParameters"]
recoParameters(mdf::MDFv2InMemory, parameters) = mdf.custom["recoParameters"] = parameters

measTemperatures(mdf::MDFv2InMemory) = @keyoptional mdf.custom["measTemperature"]
measTemperatures(mdf::MDFv2InMemory, measTemperatures) = mdf.custom["measTemperature"] = measTemperatures

measObservedDriveField(mdf::MDFv2InMemory) = @keyoptional mdf.custom["measObservedDriveField"]
measObservedDriveField(mdf::MDFv2InMemory, measDriveFields) = mdf.custom["measObservedDriveField"] = measDriveFields

measAppliedDriveField(mdf::MDFv2InMemory) = @keyoptional mdf.custom["measAppliedDriveField"]
measAppliedDriveField(mdf::MDFv2InMemory, measTransmit) = mdf.custom["measAppliedDriveField"] = measTransmit

auxiliaryData(mdf::MDFv2InMemory) = @keyoptional mdf.custom["auxiliaryData"]
auxiliaryData(mdf::MDFv2InMemory, auxiliaryData) = mdf.custom["auxiliaryData"] = auxiliaryData

export measDataRaw
function measDataRaw(mdf::MDFv2InMemory)
  if !(isnothing(mdf.measurement))
    return mdf.measurement.data
  else
    return mdf.measurement
  end
end
function measDataRaw(mdf::MDFv2InMemory, value)
  # Automatically create fields if they do not exist
  if isnothing(mdf.measurement)
    @debug "Creating field measurement"
    mdf.measurement = MDFv2Measurement()
  end
  
  mdf.measurement.data = value
end

filepath(mdf::MDFv2InMemory) = nothing # Has to be implemented...

# Helpers

export addTracer
function addTracer(mdfPart::MDFv2Tracer;
                   batch::Union{String, Missing},
                   concentration::Union{Float64, Missing},
                   injectionTime::Union{DateTime, Nothing},
                   name::Union{String, Missing},
                   solute::Union{String, Missing},
                   vendor::Union{String, Missing},
                   volume::Union{Float64, Missing})

  if !ismissing(tracerBatch(mdfPart))
    if !ismissing(batch)
      tracerBatch(mdfPart, vcat(tracerBatch(mdfPart), batch))
    else
      tracerBatch(mdfPart, vcat(tracerBatch(mdfPart), "N.A."))
    end
  else
    tracerBatch(mdfPart, [batch])
  end

  if !ismissing(tracerConcentration(mdfPart))
    if !ismissing(concentration)
      tracerConcentration(mdfPart, vcat(tracerConcentration(mdfPart), concentration))
    else
      tracerConcentration(mdfPart, vcat(tracerConcentration(mdfPart), -1.0))
    end
  else
    tracerConcentration(mdfPart, [concentration])
  end

  if !isnothing(tracerInjectionTime(mdfPart))
    if !isnothing(injectionTime)
      tracerInjectionTime(mdfPart, vcat(tracerInjectionTime(mdfPart), injectionTime))
    else
      tracerInjectionTime(mdfPart, vcat(tracerInjectionTime(mdfPart), DateTime(0)))
    end
  else
    tracerInjectionTime(mdfPart, [injectionTime])
  end

  if !ismissing(tracerName(mdfPart))
    if !ismissing(name)
      tracerName(mdfPart, vcat(tracerName(mdfPart), name))
    else
      tracerName(mdfPart, vcat(tracerName(mdfPart), "N.A."))
    end
  else
    tracerName(mdfPart, [name])
  end

  if !ismissing(tracerSolute(mdfPart))
    if !ismissing(solute)
      tracerSolute(mdfPart, vcat(tracerSolute(mdfPart), solute))
    else
      tracerSolute(mdfPart, vcat(tracerSolute(mdfPart), "N.A."))
    end
  else
    tracerSolute(mdfPart, [solute])
  end

  if !ismissing(tracerVendor(mdfPart))
    if !ismissing(vendor)
      tracerVendor(mdfPart, vcat(tracerVendor(mdfPart), vendor))
    else
      tracerVendor(mdfPart, vcat(tracerVendor(mdfPart), "N.A."))
    end
  else
    tracerVendor(mdfPart, [vendor])
  end

  if !ismissing(tracerVolume(mdfPart))
    if !ismissing(volume)
      tracerVolume(mdfPart, vcat(tracerVolume(mdfPart), volume))
    else
      tracerVolume(mdfPart, vcat(tracerVolume(mdfPart), "N.A."))
    end
  else
    tracerVolume(mdfPart, [volume])
  end

end
addTracer(mdf::MDFv2InMemory; kwargs...) = addTracer(tracer(mdf); kwargs...)