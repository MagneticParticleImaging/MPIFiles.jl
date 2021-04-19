using Unitful

export MDFv2Root, MDFv2Study, MDFv2Experiment, MDFv2Tracer, MDFv2Scanner,
       MDFv2Drivefield, MDFv2Receiver, MDFv2Acquisition,
       MDFv2Measurement, MDFv2Calibration, MDFv2Reconstruction, MDFv2InMemory

"Root group of an in-memory MDF"
mutable struct MDFv2Root
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

defaultMDFv2Root() = MDFv2Root(time=Dates.now(), uuid=UUIDs.uuid4(), "2.1.0")

"Study group of an in-memory MDF"
mutable struct MDFv2Study
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

defaultMDFv2Study() = MDFv2Study()

"Experiment group of an in-memory MDF"
mutable struct MDFv2Experiment
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

defaultMDFv2Experiment() = MDFv2Experiment() # Should we create the UUID automatically?

"Tracer group of an in-memory MDF; optional"
mutable struct MDFv2Tracer
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

"Scanner group of an in-memory MDF"
mutable struct MDFv2Scanner
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

"Drivefield subgroup of acquisition group of an in-memory MDF"
mutable struct MDFv2Drivefield
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

"Receiver subgroup of acquisition group of an in-memory MDF"
mutable struct MDFv2Receiver
  "Bandwidth of the receiver unit"
  bandwidth::Union{Float64, Missing}
  "Dimension less scaling factor and offset (a_c, b_c) to convert raw data into a
  physical quantity with corresponding unit of measurement unit; optional"
  dataConversionFactor::Union{Array{Float64, 2}, Nothing}
  "Induction factor mapping the projection of the magnetic moment to the voltage in the receive coil; optional"
  inductionFactor::Union{Array{Float64, 2}, Nothing}
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
    dataConversionFactor = missing,
    inductionFactor = missing,
    numChannels = missing,
    numSamplingPoints = missing,
    transferFunction = missing,
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

"Acquisition group of an in-memory MDF"
mutable struct MDFv2Acquisition
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
  time::Union{DateTime, Missing}

  drivefield::MDFv2Drivefield
  receiver::MDFv2Receiver

  function MDFv2Acquisition(;
    gradient = nothing,
    numAverages = missing,
    numFrames = missing,
    numPeriodsPerFrame = missing,
    offsetField = missing,
    time = missing,
    drivefield = nothing,
    receiver = missing)
    
    return new(
      gradient,
      numAverages,
      numFrames,
      numPeriodsPerFrame,
      offsetField,
      time,
      drivefield,
      receiver
    )
  end
end

defaultMDFv2Acquisition() = MDFv2Acquisition()

"Measurement group of an in-memory MDF"
mutable struct MDFv2Measurement
  "Measured data at a specific processing stage"
  data::Union{Array{Number, 4}, Missing} # Should be mmapable
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

"Calibration group of an in-memory MDF"
mutable struct MDFv2Calibration
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
  size::Union{Vector{Float64}, Nothing}
  "Signal-to-noise estimate for recorded frequency components; optional"
  snr::Union{Array{Float64, 3}, Nothing}

  function MDFv2Calibration(;
    deltaSampleSize = nothing,
    fieldOfView = nothing,
    fieldOfViewCenter = nothing,
    method = missing,
    offsetFields = nothing,
    order = nothing,
    positions = nothing,
    size = nothing,
    snr = nothing)

    return new(
      deltaSampleSize,
      fieldOfView,
      fieldOfViewCenter,
      method,
      offsetFields,
      order,
      positions,
      size,
      snr
    )
  end
end

defaultMDFv2Calibration() = MDFv2Calibration(order="xyz")

"Reconstruction group of an in-memory MDF"
mutable struct MDFv2Reconstruction
  "Reconstructed data"
  data::Union{Array{Number, 3}, Missing}
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

  function MDFv2InMemory(;
    root = missing,
    study = missing,
    experiment = missing,
    tracer = nothing,
    scanner = missing,
    acquisition = missing,
    measurement = nothing,
    calibration = nothing,
    reconstruction = nothing)

    return new(
      root,
      study,
      experiment,
      tracer,
      scanner,
      acquisition,
      measurement,
      calibration,
      reconstruction
    )
  end
end

function defaultMDFv2InMemory()
  return MDFv2InMemory(
    root = defaultMDFv2Root(),
    study = defaultMDFv2Study(),
    experiment = defaultMDFv2Experiment(),
    scanner = defaultMDFv2Scanner(),
    acquisition = defaultMDFv2Acquisition()
  )
end

"Check, whether all non-optional fields have been
set and if the dimensions of the fields match"
function checkConsistency(mdf::MDFv2InMemory)
  for field in fieldnames(MDFv2InMemory)
    
  end
end

  
# # general parameters
# version(f::MDFFile)::VersionNumber = VersionNumber( f["/version"] )
# uuid(f::MDFFile)::UUID = UUID(f["/uuid"])
# time(f::MDFFileV1)::DateTime = DateTime( f["/date"] )
# time(f::MDFFileV2)::DateTime = DateTime( f["/time"] )

# # study parameters
# studyName(f::MDFFile)::String = f["/study/name"]
# studyNumber(f::MDFFileV1)::Int = 0
# studyNumber(f::MDFFileV2)::Int = f["/study/number"]
# studyUuid(f::MDFFileV1) = nothing
# studyUuid(f::MDFFileV2) = UUID(f["/study/uuid"])
# studyDescription(f::MDFFileV1)::String = "n.a."
# studyDescription(f::MDFFileV2)::String = f["/study/description"]
# function studyTime(f::MDFFile)
#   t = f["/study/time"]
#   if typeof(t)==String
#     return DateTime(t)
#   else
#     return nothing
#   end
# end

# # experiment parameters
# experimentName(f::MDFFileV1)::String = "n.a."
# experimentName(f::MDFFileV2)::String = f["/experiment/name"]
# experimentNumber(f::MDFFileV1)::Int64 = parse(Int64, f["/study/experiment"])
# experimentNumber(f::MDFFileV2)::Int64 = f["/experiment/number"]
# experimentUuid(f::MDFFileV1) = nothing
# experimentUuid(f::MDFFileV2) = UUID(f["/experiment/uuid"])
# experimentDescription(f::MDFFileV1)::String = f["/study/description"]
# experimentDescription(f::MDFFileV2)::String = f["/experiment/description"]
# experimentSubject(f::MDFFileV1)::String = f["/study/subject"]
# experimentSubject(f::MDFFileV2)::String = f["/experiment/subject"]
# experimentIsSimulation(f::MDFFileV2)::Bool = Bool( f["/experiment/isSimulation"] )
# experimentIsSimulation(f::MDFFileV1)::Bool = Bool( f["/study/simulation"] )
# experimentIsCalibration(f::MDFFile)::Bool = haskey(f.file, "/calibration")
# experimentHasReconstruction(f::MDFFile)::Bool = haskey(f.file, "/reconstruction")
# experimentHasMeasurement(f::MDFFileV1)::Bool = haskey(f.file, "/measurement") ||
#                                           haskey(f.file, "/calibration")
# experimentHasMeasurement(f::MDFFileV2)::Bool = haskey(f.file, "/measurement")

# _makeStringArray(s::String) = [s]
# _makeStringArray(s::Vector{T}) where {T<:AbstractString} = s

# # tracer parameters
# tracerName(f::MDFFileV1)::Vector{String} = [f["/tracer/name"]]
# tracerName(f::MDFFileV2)::Vector{String} = _makeStringArray(f["/tracer/name"])
# tracerBatch(f::MDFFileV1)::Vector{String} = [f["/tracer/batch"]]
# tracerBatch(f::MDFFileV2)::Vector{String} = _makeStringArray(f["/tracer/batch"])
# tracerVolume(f::MDFFileV1)::Vector{Float64} = [f["/tracer/volume"]]
# tracerVolume(f::MDFFileV2)::Vector{Float64} = [f["/tracer/volume"]...]
# tracerConcentration(f::MDFFileV1)::Vector{Float64} = [f["/tracer/concentration"]]
# tracerConcentration(f::MDFFileV2)::Vector{Float64} = [f["/tracer/concentration"]...]
# tracerSolute(f::MDFFileV2)::Vector{String} = _makeStringArray(f["/tracer/solute"])
# tracerSolute(f::MDFFileV1)::Vector{String} = ["Fe"]
# function tracerInjectionTime(f::MDFFile)::Vector{DateTime}
#   p = typeof(f) <: MDFFileV1 ? "/tracer/time" : "/tracer/injectionTime"
#   if f[p] == nothing
#     return nothing
#   end

#   if typeof(f[p]) == String
#     return [DateTime(f[p])]
#   else
#     return [DateTime(y) for y in f[p]]
#   end
# end
# #tracerInjectionTime(f::MDFFileV2) = DateTime( f["/tracer/injectionTime"] )
# tracerVendor(f::MDFFileV1)::Vector{String} = [f["/tracer/vendor"]]
# tracerVendor(f::MDFFileV2)::Vector{String} = _makeStringArray(f["/tracer/vendor"])

# # scanner parameters
# scannerFacility(f::MDFFile)::String = f["/scanner/facility"]
# scannerOperator(f::MDFFile)::String = f["/scanner/operator"]
# scannerManufacturer(f::MDFFile)::String = f["/scanner/manufacturer"]
# scannerName(f::MDFFileV1)::String = f["/scanner/model"]
# scannerName(f::MDFFileV2)::String = f["/scanner/name", ""]
# scannerTopology(f::MDFFile)::String = f["/scanner/topology"]

# # acquisition parameters
# acqStartTime(f::MDFFileV1)::DateTime = DateTime( f["/acquisition/time"] )
# acqStartTime(f::MDFFileV2)::DateTime = DateTime( f["/acquisition/startTime"] )
# acqNumAverages(f::MDFFileV1)::Int = f["/acquisition/drivefield/averages"]
# acqNumAverages(f::MDFFileV2)::Int = f["/acquisition/numAverages",1]
# function acqNumFrames(f::MDFFileV1)::Int
#   if experimentIsCalibration(f)
#     return size(f.mmap_measData,2)
#   else
#     return f["/acquisition/numFrames"]
#   end
# end
# acqNumFrames(f::MDFFileV2)::Int = f["/acquisition/numFrames"]
# acqNumPeriodsPerFrame(f::MDFFileV1)::Int = 1
# acqNumPeriodsPerFrame(f::MDFFileV2)::Int = f["/acquisition/numPeriods",1]

# acqGradient(f::MDFFileV1)::Array{Float64,4} = reshape(Matrix(Diagonal(f["/acquisition/gradient"])), 3,3,1,1)
# function acqGradient(f::MDFFileV2)::Array{Float64,4} 
#   G = f["/acquisition/gradient"]
#   if ndims(G) == 4
#     return G
#   elseif ndims(G) == 3 # for corrupt files
#     return reshape(G,3,3,1,size(G,3))
#   elseif ndims(G) == 2 && prod(size(G)) == 9  # for corrupt files
#     return reshape(G,3,3,1,1)
#   else # for corrupt files
#     return reshape(Matrix(Diagonal(vec(G))),3,3,1,1)
#   end
# end

# acqOffsetField(f::MDFFileV1)::Array{Float64,3} = f["/acquisition/offsetField", reshape([0.0,0.0,0.0],3,1,1)  ]
# function acqOffsetField(f::MDFFileV2)::Array{Float64,3} 
#   H = f["/acquisition/offsetField", reshape([0.0,0.0,0.0],3,1,1)  ]
#   if ndims(H) == 3
#     return H
#   else # for corrupt files
#     return reshape(H,:,1,1)
#   end
# end

# # drive-field parameters
# dfNumChannels(f::MDFFile)::Int = f["/acquisition/drivefield/numChannels"]
# dfStrength(f::MDFFileV1)::Array{Float64,3} = addTrailingSingleton( addLeadingSingleton(
#           f["/acquisition/drivefield/strength"], 2), 3)
# dfStrength(f::MDFFileV2)::Array{Float64,3} = f["/acquisition/drivefield/strength"]
# dfPhase(f::MDFFileV1)::Array{Float64,3} = dfStrength(f) .*0 .+  1.5707963267948966 # Bruker specific!
# dfPhase(f::MDFFileV2)::Array{Float64,3} = f["/acquisition/drivefield/phase"]
# dfBaseFrequency(f::MDFFile)::Float64 = f["/acquisition/drivefield/baseFrequency"]
# dfCustomWaveform(f::MDFFileV2)::String = f["/acquisition/drivefield/customWaveform"]
# dfDivider(f::MDFFileV1) = addTrailingSingleton(
#                 f["/acquisition/drivefield/divider"],2)
# dfDivider(f::MDFFileV2) = f["/acquisition/drivefield/divider"]
# dfWaveform(f::MDFFileV1)::String = "sine"
# dfWaveform(f::MDFFileV2)::String = f["/acquisition/drivefield/waveform"]
# dfCycle(f::MDFFile)::Float64 = f["/acquisition/drivefield/cycle"]
# dfCycle(f::MDFFileV1)::Float64 = f["/acquisition/drivefield/period"]

# # receiver parameters
# rxNumChannels(f::MDFFile)::Int64 = f["/acquisition/receiver/numChannels"]
# rxBandwidth(f::MDFFile)::Float64 = f["/acquisition/receiver/bandwidth"]
# rxNumSamplingPoints(f::MDFFile)::Int64 = f["/acquisition/receiver/numSamplingPoints"]
# function rxTransferFunction(f::MDFFile)
#   parameter = "/acquisition/receiver/transferFunction"
#   if haskey(f.file, parameter)
#     return readComplexArray(f.filename, parameter)
#   else
#     return nothing
#   end
# end
# function rxTransferFunctionFileName(f::MDFFile)
#   parameter = "/acquisition/receiver/transferFunctionFileName"
#   if haskey(f.file, parameter)
#     return f[parameter]
#   else
#     return nothing
#   end
# end
# function rxHasTransferFunction(f::MDFFile)
#   haskey(f.file, "/acquisition/receiver/transferFunction")
# end
# rxInductionFactor(f::MDFFileV1) = nothing
# rxInductionFactor(f::MDFFileV2) = f["/acquisition/receiver/inductionFactor"]

# rxUnit(f::MDFFileV1)::String = "a.u."
# rxUnit(f::MDFFileV2)::String = f["/acquisition/receiver/unit"]
# rxDataConversionFactor(f::MDFFileV1) = repeat([1.0, 0.0], outer=(1,rxNumChannels(f)))
# rxDataConversionFactor(f::MDFFileV2) = f["/acquisition/receiver/dataConversionFactor"]

# # measurements
# function measData(f::MDFFileV1, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
#                   receivers=1:rxNumChannels(f))
#   if !haskey(f.file, "/measurement")
#     # the V1 file is a calibration
#     data = f["/calibration/dataFD"]
#     if ndims(data) == 4
#       return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4),1))
#     else
#       return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4),size(data,5)))
#     end
#   end
#   tdExists = haskey(f.file, "/measurement/dataTD")

#   if tdExists
#     data = zeros(Float64, rxNumSamplingPoints(f), length(receivers), length(frames))
#     for (i,fr) in enumerate(frames)
#       data[:,:,:,i] = f.mmap_measData[:, receivers, fr]
#     end
#     return reshape(data,size(data,1),size(data,2),1,size(data,3))
#   else
#     data = zeros(Float64, 2, rxNumFrequencies(f), length(receivers), length(frames))
#     for (i,fr) in enumerate(frames)
#       data[:,:,:,i] = f.mmap_measData[:,:,receivers, fr]
#     end

#     dataFD = reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4)))
#     dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
#     return reshape(dataTD,size(dataTD,1),size(dataTD,2),1,size(dataTD,3))
#   end
# end

# function measData(f::MDFFileV2, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
#                   receivers=1:rxNumChannels(f))

#   if measIsFastFrameAxis(f)
#     data = f.mmap_measData[frames, :, receivers, periods]
#     data = reshape(data, length(frames), size(f.mmap_measData,2), length(receivers), length(periods))
#   else
#     data = f.mmap_measData[:, receivers, periods, frames]
#     data = reshape(data, size(f.mmap_measData,1), length(receivers), length(periods), length(frames))
#   end
#   return data
# end


# function measDataTDPeriods(f::MDFFileV1, periods=1:acqNumPeriods(f),
#                   receivers=1:rxNumChannels(f))
#   tdExists = haskey(f.file, "/measurement/dataTD")

#   if tdExists
#     data = f.mmap_measData[:, receivers, periods]
#     return data
#   else
#     data = f.mmap_measData[:, :, receivers, periods]

#     dataFD = reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4)))
#     dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
#     return dataTD
#   end
# end


# function measDataTDPeriods(f::MDFFileV2, periods=1:acqNumPeriods(f),
#                   receivers=1:rxNumChannels(f))
#   if measIsFastFrameAxis(f)
#     error("measDataTDPeriods can currently not handle transposed data!")
#   end

#   data = reshape(f.mmap_measData,Val(3))[:, receivers, periods]

#   return data
# end

# function systemMatrix(f::MDFFileV1, rows, bgCorrection=true)
#   if !experimentIsCalibration(f)
#     return nothing
#   end

#   data = reshape(f.mmap_measData,Val(3))[:, :, rows]
#   return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3)))
# end

# function systemMatrix(f::MDFFileV2, rows, bgCorrection=true)
#   if !haskey(f.file, "/measurement") || !measIsFastFrameAxis(f) ||
#     !measIsFourierTransformed(f)
#     return nothing
#   end

#   rows_ = rowsToSubsampledRows(f, rows)

#   data_ = reshape(f.mmap_measData, size(f.mmap_measData,1),
#                                     size(f.mmap_measData,2)*size(f.mmap_measData,3),
#                                     size(f.mmap_measData,4))[:, rows_, :]
#   data = reshape(data_, Val(2))

#   fgdata = data[measFGFrameIdx(f),:]

#   if measIsSparsityTransformed(f)
#     dataBackTrafo = similar(fgdata, prod(calibSize(f)), size(fgdata,2))
#     B = linearOperator(f["/measurement/sparsityTransformation"], calibSize(f))

#     tmp = f["/measurement/subsamplingIndices"]
#     subsamplingIndices_ = reshape(tmp, size(tmp,1),
#                                       size(tmp,2)*size(tmp,3),
#                                       size(tmp,4))[:, rows_, :]
#     subsamplingIndices = reshape(subsamplingIndices_, Val(2))

#     for l=1:size(fgdata,2)
#       dataBackTrafo[:,l] .= 0.0
#       dataBackTrafo[subsamplingIndices[:,l],l] .= fgdata[:,l]
#       dataBackTrafo[:,l] .= adjoint(B) * vec(dataBackTrafo[:,l])
#     end
#     fgdata = dataBackTrafo
#   end

#   if bgCorrection # this assumes equidistent bg frames
#     @debug "Applying bg correction on system matrix (MDF)"
#     bgdata = data[measBGFrameIdx(f),:]
#     blockLen = measBGFrameBlockLengths( invpermute!(measIsBGFrame(f), measFramePermutation(f)) )
#     st = 1
#     for j=1:length(blockLen)
#       bgdata[st:st+blockLen[j]-1,:] .=
#             mean(bgdata[st:st+blockLen[j]-1,:], dims=1)
#       st += blockLen[j]
#     end

#     bgdataInterp = interpolate(bgdata, (BSpline(Linear()), NoInterp()))
#     # Cubic does not work for complex numbers
#     origIndex = measFramePermutation(f)
#     M = size(fgdata,1)
#     K = size(bgdata,1)
#     N = M + K
#     for m=1:M
#       alpha = (origIndex[m]-1)/(N-1)*(K-1)+1
#       for k=1:size(fgdata,2)
#         fgdata[m,k] -= bgdataInterp(alpha,k)
#       end
#     end
#   end
#   return fgdata
# end

# function systemMatrixWithBG(f::MDFFileV2)
#   if !haskey(f.file, "/measurement") || !measIsFastFrameAxis(f) ||
#       !measIsFourierTransformed(f)
#       return nothing
#   end

#   data = f.mmap_measData[:, :, :, :]
#   return data
# end

# # This is a special variant used for matrix compression
# function systemMatrixWithBG(f::MDFFileV2, freq)
#   if !haskey(f.file, "/measurement") || !measIsFastFrameAxis(f) ||
#     !measIsFourierTransformed(f)
#     return nothing
#   end

#   data = f.mmap_measData[:, freq, :, :]
#   return data
# end

# function measIsFourierTransformed(f::MDFFileV1)
#   if !experimentIsCalibration(f)
#     return false
#   else
#     return true
#   end
# end
# measIsFourierTransformed(f::MDFFileV2) = Bool(f["/measurement/isFourierTransformed"])

# measIsTFCorrected(f::MDFFileV1) = false
# measIsTFCorrected(f::MDFFileV2) = Bool(f["/measurement/isTransferFunctionCorrected"])

# measIsSpectralLeakageCorrected(f::MDFFileV1) = false
# measIsSpectralLeakageCorrected(f::MDFFileV2) = Bool(f["/measurement/isSpectralLeakageCorrected"])

# function measIsBGCorrected(f::MDFFileV1)
#   if !experimentIsCalibration(f)
#     return false
#   else
#     return true
#   end
# end
# measIsBGCorrected(f::MDFFileV2) = Bool(f["/measurement/isBackgroundCorrected"])

# measIsFrequencySelection(f::MDFFileV1) = false
# measIsFrequencySelection(f::MDFFileV2) = Bool(f["/measurement/isFrequencySelection"])
# measFrequencySelection(f::MDFFileV2) = f["/measurement/frequencySelection"]

# measIsSparsityTransformed(f::MDFFileV1) = false
# function measIsSparsityTransformed(f::MDFFileV2)
#   if haskey(f.file, "/measurement/isSparsityTransformed")
#     Bool(f["/measurement/isSparsityTransformed"])
#   else
#     return false
#   end
# end

# function measIsFastFrameAxis(f::MDFFileV1)
#   if !experimentIsCalibration(f)
#     return false
#   else
#     return true
#   end
# end

# function measIsFastFrameAxis(f::MDFFileV2)
#   if haskey(f.file, "/measurement/isFastFrameAxis")
#     return Bool(f["/measurement/isFastFrameAxis"])
#   else
#     @warn "/measurement/isFastFrameAxis missing in MDF data set. `measIsFastFrameAxis` returning false per default."
#     return false
#   end
# end

# function measIsFramePermutation(f::MDFFileV1)
#   if !experimentIsCalibration(f)
#     return false
#   else
#     return true
#   end
# end
# measIsFramePermutation(f::MDFFileV2) = Bool(f["/measurement/isFramePermutation"])
# measIsBGFrame(f::MDFFileV1) = zeros(Bool, acqNumFrames(f))
# measIsBGFrame(f::MDFFileV2) = convert(Array{Bool},f["/measurement/isBackgroundFrame"])
# measFramePermutation(f::MDFFileV1) = nothing
# measFramePermutation(f::MDFFileV2) = f["/measurement/framePermutation"]
# fullFramePermutation(f::MDFFile) = fullFramePermutation(f, calibIsMeanderingGrid(f))

# measIsCalibProcessed(f::MDFFile) = measIsFramePermutation(f) && 
#                                     measIsFourierTransformed(f) &&
#                                     measIsFastFrameAxis(f)

# #calibrations
# calibSNR(f::MDFFileV1) = addTrailingSingleton(f["/calibration/snrFD"],3)
# calibSNR(f::MDFFileV2) = f["/calibration/snr"]
# calibFov(f::MDFFile) = f["/calibration/fieldOfView"]
# calibFovCenter(f::MDFFile) = f["/calibration/fieldOfViewCenter"]
# calibSize(f::MDFFile) = f["/calibration/size"]
# calibOrder(f::MDFFile) = f["/calibration/order"]
# calibOffsetField(f::MDFFile) = f["/calibration/offsetField"]
# calibDeltaSampleSize(f::MDFFile) = f["/calibration/deltaSampleSize",[0.0,0.0,0.0]]
# calibMethod(f::MDFFile) = f["/calibration/method"]
# calibIsMeanderingGrid(f::MDFFile) = Bool(f["/calibration/isMeanderingGrid", 0])
# calibPositions(f::MDFFile) = f["/calibration/positions"]

# # reconstruction results
# recoData(f::MDFFileV1) = addLeadingSingleton(
#           f[ "/reconstruction/data"], 3)
# recoData(f::MDFFileV2) = f["/reconstruction/data"]
# recoFov(f::MDFFile)::Vector{Float64} = f["/reconstruction/fieldOfView"]
# recoFovCenter(f::MDFFile)::Vector{Float64} = f["/reconstruction/fieldOfViewCenter"]
# recoSize(f::MDFFile)::Vector{Int64} = f["/reconstruction/size"]
# recoOrder(f::MDFFile) = f["/reconstruction/order"]
# recoPositions(f::MDFFile) = f["/reconstruction/positions"]

# # this is non-standard
# function recoParameters(f::MDFFile)
#   if !haskey(f.file, "/reconstruction/_parameters")
#     return nothing
#   else
#     return loadParams(f.file, "/reconstruction/_parameters")
#   end
# end

# # additional functions that should be implemented by an MPIFile
# filepath(f::MDFFile) = "RAM" # Yes, this is a bit weird, but has to be implemented