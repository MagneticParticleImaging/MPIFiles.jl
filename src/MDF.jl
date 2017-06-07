using HDF5

export MDFFile, MDFFileV1, MDFFileV2

abstract MDFFile <: MPIFile

# We use a dedicated type for v1 and v2. If both versions
# are the same we use the abstract type MDFFile
type MDFFileV1 <: MDFFile
  filename::String
  param_cache
  mmap_measData
  mmap_procData
end

MDFFileV1(filename::String) = MDFFileV1(filename,Dict{String,Any}(),nothing,nothing)

type MDFFileV2 <: MDFFile
  filename::String
  param_cache
  mmap_measData
  mmap_procData
end

MDFFileV2(filename::String) = MDFFileV2(filename,Dict{String,Any}(),nothing,nothing)

# This dispatches on the file extension and automatically
# generates the correct type
function (::Type{MDFFile})(filename::String)
  vers = VersionNumber( h5read(filename, "/version") )
  if vers < v"2.0"
    return MDFFileV1(filename)
  else
    return MDFFileV2(filename)
  end
end

function Base.show(io::IO, f::MDFFileV1)
  print(io, "MDF v1: ", f.filename)
end

function Base.show(io::IO, f::MDFFileV2)
  print(io, "MDF v2: ", f.filename)
end

function h5exists(filename, parameter)
  return h5open(filename) do file
    exists(file, parameter)
  end
end

function h5readornull(filename, parameter)
  if h5exists(filename, parameter)
    return h5read(filename, parameter)
  else
    return nothing
  end
end

function getindex(f::MDFFile, parameter)
  if !haskey(f.param_cache,parameter)
    f.param_cache[parameter] = h5readornull(f.filename, parameter)
  end
  return f.param_cache[parameter]
end



# general parameters
version(f::MDFFile) = VersionNumber( f["/version"] )
uuid(f::MDFFile) = f["/uuid"]
time(f::MDFFileV1) = DateTime( f["/date"] )
time(f::MDFFileV2) = DateTime( f["/time"] )

# study parameters
studyName(f::MDFFile) = f["/study/name"]
studyNumber(f::MDFFileV1) = 0
studyNumber(f::MDFFileV2) = f["/study/number"]
studyDescription(f::MDFFileV1) = "n.a."
studyDescription(f::MDFFileV2) = f["/study/description"]

# experiment parameters
experimentName(f::MDFFileV1) = "n.a."
experimentName(f::MDFFileV2) = f["/experiment/name"]
experimentNumber(f::MDFFileV1) = parse(Int64, f["/study/experiment"])
experimentNumber(f::MDFFileV2) = f["/experiment/number"]
experimentDescription(f::MDFFileV1) = f["/study/description"]
experimentDescription(f::MDFFileV2) = f["/experiment/description"]
experimentSubject(f::MDFFileV1) = f["/study/subject"]
experimentSubject(f::MDFFileV2) = f["/experiment/subject"]
experimentIsSimulation(f::MDFFileV2) = Bool( f["/experiment/isSimulation"] )
experimentIsSimulation(f::MDFFileV1) = Bool( f["/study/simulation"] )
experimentIsCalibration(f::MDFFile) = h5exists(f.filename, "/calibration")
experimentHasProcessing(f::MDFFileV1) = experimentIsCalibration(f)
experimentHasProcessing(f::MDFFileV2) = h5exists(f.filename, "/processing")
experimentHasReconstruction(f::MDFFile) = h5exists(f.filename, "/reconstruction")
experimentHasMeasurement(f::MDFFile) = h5exists(f.filename, "/measurement")
# tracer parameters
tracerName(f::MDFFile) = f["/tracer/name"]
tracerBatch(f::MDFFile) = f["/tracer/batch"]
tracerVolume(f::MDFFile) = f["/tracer/volume"]
tracerConcentration(f::MDFFile) = f["/tracer/concentration"]
tracerSolute(f::MDFFileV2) = f["/tracer/solute"]
tracerSolute(f::MDFFileV1) = "Fe"
tracerInjectionTime(f::MDFFileV1) = DateTime( f["/tracer/time"] )
tracerInjectionTime(f::MDFFileV2) = DateTime( f["/tracer/injectionTime"] )
tracerVendor(f::MDFFile) = f["/tracer/vendor"]

# scanner parameters
scannerFacility(f::MDFFile) = f["/scanner/facility"]
scannerOperator(f::MDFFile) = f["/scanner/operator"]
scannerManufacturer(f::MDFFile) = f["/scanner/manufacturer"]
scannerModel(f::MDFFile) = f["/scanner/model"]
scannerTopology(f::MDFFile) = f["/scanner/topology"]

# acquisition parameters
acqStartTime(f::MDFFileV1) = DateTime( f["/acquisition/time"] )
acqStartTime(f::MDFFileV2) = DateTime( f["/acquisition/startTime"] )
acqNumFrames(f::MDFFile) = f["/acquisition/numFrames"]
acqNumBGFrames(f::MDFFileV1) = 0
acqNumBGFrames(f::MDFFileV2) = f["/acquisition/numBackgroundFrames"]
acqFramePeriod(f::MDFFile) = f["/acquisition/framePeriod"]
acqNumPatches(f::MDFFile) = f["/acquisition/numPatches"]
acqGradient(f::MDFFileV1) = addLeadingSingleton(f["/acquisition/gradient"],2)
acqGradient(f::MDFFileV2) = f["/acquisition/gradient"]
acqOffsetField(f::MDFFile) = f["/acquisition/offsetField"]
acqOffsetFieldShift(f::MDFFileV1) = addLeadingSingleton(
              f["/acquisition/drivefield/fieldOfViewCenter"],2 )
acqOffsetFieldShift(f::MDFFileV2) = f["/acquisition/offsetFieldShift"]

# drive-field parameters
dfNumChannels(f::MDFFile) = f["/acquisition/drivefield/numChannels"]
dfStrength(f::MDFFileV1) = addTrailingSingleton( addLeadingSingleton(
         f["/acquisition/drivefield/strength"], 2), 3)
dfStrength(f::MDFFileV2) = f["/acquisition/drivefield/strength"]
dfPhase(f::MDFFileV1) = dfStrength(f) .*0 .+  1.5707963267948966 # Bruker specific!
dfPhase(f::MDFFileV2) = f["/acquisition/drivefield/phase"]
dfBaseFrequency(f::MDFFile) = f["/acquisition/drivefield/baseFrequency"]
dfCustomWaveform(f::MDFFileV2) = f["/acquisition/drivefield/customWaveform"]
dfDivider(f::MDFFileV1) = addTrailingSingleton(
                f["/acquisition/drivefield/divider"],2)
dfDivider(f::MDFFileV2) = f["/acquisition/drivefield/divider"]
dfWaveform(f::MDFFileV1) = "sine"
dfWaveform(f::MDFFileV2) = f["/acquisition/drivefield/waveform"]
dfPeriod(f::MDFFile) = f["/acquisition/drivefield/period"]

# receiver parameters
rxNumChannels(f::MDFFile) = f["/acquisition/receiver/numChannels"]
rxNumAverages(f::MDFFileV1) = f["/acquisition/drivefield/averages"]
rxNumAverages(f::MDFFileV2) = f["/acquisition/receiver/numAverages"]
rxBandwidth(f::MDFFile) = f["/acquisition/receiver/bandwidth"]
rxNumSamplingPoints(f::MDFFile) = f["/acquisition/receiver/numSamplingPoints"]
rxTransferFunction(f::MDFFile) = f["/acquisition/receiver/transferFunction"]

# measurements
measUnit(f::MDFFileV1) = "a.u."
measUnit(f::MDFFileV2) = f["/measurement/unit"]
measDataConversionFactor(f::MDFFileV1) = [1.0, 0.0]
measDataConversionFactor(f::MDFFileV2) = f["/measurement/dataConversionFactor"]
function measData(f::MDFFileV1, frames=1:acqNumAllFrames(f), patches=1:acqNumPatches(f),
                  receivers=1:rxNumChannels(f))
  if !h5exists(f.filename, "/measurement")
    return nothing
  end
  tdExists = h5exists(f.filename, "/measurement/dataTD")

  if tdExists
    if f.mmap_measData == nothing
      h5open(f.filename,"r") do file
        f.mmap_measData = readmmap(file["/measurement/dataTD"])
      end
    end
    #data = zeros(Float64, dims[1], length(receivers), length(frames))
    data = zeros(Float64, rxNumSamplingPoints(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:, receivers, fr]
      #h5read(f.filename, "/measurement/dataTD", (:,  receivers, fr) )
    end
    return reshape(data,size(data,1),size(data,2),1,size(data,3))
  else
    if f.mmap_measData == nothing
      h5open(f.filename,"r") do file
        f.mmap_measData = readmmap(file["/measurement/dataFD"])
      end
    end
    #data = zeros(Float64, dims[1], dims[2], length(receivers), length(frames))
    data = zeros(Float64, 2, rxNumFrequencies(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:,:,receivers, fr]
      #h5read(f.filename, "/measurement/dataFD", (:, :, receivers, fr) )
    end

    dataFD = reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4)))
    dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
    return reshape(dataTD,size(dataTD,1),size(dataTD,2),1,size(dataTD,3))
  end
end
function measData(f::MDFFileV2, frames=1:acqNumAllFrames(f), patches=1:acqNumPatches(f),
                  receivers=1:rxNumChannels(f))
  if !h5exists(f.filename, "/measurement")
    return nothing
  end
  if f.mmap_measData == nothing
    h5open(f.filename,"r") do file
      f.mmap_measData = readmmap(file["/measurement/data"])
    end
  end
  data = zeros(Float64, rxNumSamplingPoints(f), length(receivers),
                        length(patches), length(frames))
  for (i,fr) in enumerate(frames)
    data[:,:,:,i] = f.mmap_measData[:, receivers, patches, fr]
    #h5read(f.filename, "/measurement/data", (:, receivers, patches, fr) )
  end
  return data
end
measIsBG(f::MDFFileV1) = zeros(Bool, acqNumFrames(f))
measIsBG(f::MDFFileV2) = convert(Array{Bool},f["/measurement/isBackgroundData"])

# processings
function procData(f::MDFFileV1; frames=nothing)
  if !experimentIsCalibration(f)
    return nothing
  end

  data = f["/calibration/dataFD"]
  if ndims(data) == 4
    return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),1))
  else
    return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),size(data,5)))
  end
end

function procData(f::MDFFileV1, rows)
  if !experimentIsCalibration(f)
    return nothing
  end
  if f.mmap_procData == nothing
    h5open(f.filename,"r") do file
      f.mmap_procData = readmmap(file["/calibration/dataFD"])
    end
  end
  data = f.mmap_procData[:, :, rows]
  return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3)))
end

function procData(f::MDFFileV2; frames=:)
  if !h5exists(f.filename, "/processing")
    return nothing
  end
  if f.mmap_procData == nothing
    h5open(f.filename,"r") do file
      f.mmap_procData = readmmap(file["/processing/data"])
    end
  end
  data = f.mmap_procData[:, :, :, :, frames]

  if procIsFourierTransformed(f)
    if procIsTransposed(f)
      data = f.mmap_procData[:, frames, :, :, :]
    else
      data = f.mmap_procData[:, :, :, :, frames]
    end

    return reinterpret(Complex{eltype(data)}, data,
               (size(data,2),size(data,3),size(data,4),size(data,5)))
  else
    if procIsTransposed(f)
      data = f.mmap_procData[frames, :, :, :]
    else
      data = f.mmap_procData[:, :, :, frames]
    end
    return data
  end
end

function procData(f::MDFFileV2, rows)
  if !h5exists(f.filename, "/processing") || !procIsTransposed(f)
    return nothing
  end
  if f.mmap_procData == nothing
    h5open(f.filename,"r") do file
      f.mmap_procData = readmmap(file["/processing/data"])
    end
  end
  if procIsFourierTransformed(f)
    data = f.mmap_procData[:, :, rows]
    return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3)))
  else
    data = f.mmap_procData[:, rows]
    return data
  end
end


function procIsFourierTransformed(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return true
  end
end
procIsFourierTransformed(f::MDFFileV2) = Bool(f["/processing/isFourierTransformed"])

function procIsTFCorrected(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return false
  end
end
procIsTFCorrected(f::MDFFileV2) = Bool(f["/processing/isTransferFunctionCorrected"])

function procIsAveraged(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return false
  end
end
procIsAveraged(f::MDFFileV2) = Bool(f["/processing/isAveraged"])

function procIsFramesSelected(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return false
  end
end
procIsFramesSelected(f::MDFFileV2) = Bool(f["/processing/isFramesSelected"])

function procIsBGCorrected(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return true
  end
end
procIsBGCorrected(f::MDFFileV2) = Bool(f["/processing/isBackgroundCorrected"])

function procIsTransposed(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return true
  end
end
procIsTransposed(f::MDFFileV2) = Bool(f["/processing/isTransposed"])

function procFramePermutation(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return nothing
  else
    return nothing # TODO
  end
end
procFramePermutation(f::MDFFileV2) = f["/processing/framePermutation"]


#calibrations
calibSNR(f::MDFFileV1) = addTrailingSingleton(f["/calibration/snrFD"],3)
calibSNR(f::MDFFileV2) = f["/calibration/snr"]
calibFov(f::MDFFile) = f["/calibration/fieldOfView"]
calibFovCenter(f::MDFFile) = f["/calibration/fieldOfViewCenter"]
calibSize(f::MDFFile) = f["/calibration/size"]
calibOrder(f::MDFFile) = f["/calibration/order"]
calibPositions(f::MDFFile) = f["/calibration/positions"]
calibOffsetField(f::MDFFile) = f["/calibration/offsetField"]
calibDeltaSampleSize(f::MDFFile) = f["/calibration/deltaSampleSize"]
calibMethod(f::MDFFile) = f["/calibration/method"]

# reconstruction results
recoData(f::MDFFileV1) = addLeadingSingleton(
         f[ "/reconstruction/data"], 3)
recoData(f::MDFFileV2) = f["/reconstruction/data"]
recoFov(f::MDFFile) = f["/reconstruction/fieldOfView"]
recoFovCenter(f::MDFFile) = f["/reconstruction/fieldOfViewCenter"]
recoSize(f::MDFFile) = f["/reconstruction/size"]
recoOrder(f::MDFFile) = f["/reconstruction/order"]
recoPositions(f::MDFFile) = f["/reconstruction/positions"]

# this is non-standard
function recoParameters(f::MDFFile)
  if !h5exists(f.filename, "/reconstruction/parameters")
    return nothing
  end
  return loadParams(f.filename, "/reconstruction/parameters")
end

# additional functions that should be implemented by an MPIFile
filepath(f::MDFFile) = f.filename


# Helper functions
function addLeadingSingleton(a::Array,dim)
  if ndims(a) == dim
    return a
  else
    return reshape(a,1,size(a)...)
  end
end

function addTrailingSingleton(a::Array,dim)
  if ndims(a) == dim
    return a
  else
    return reshape(a,size(a)...,1)
  end
end












#= TODO
type MDFTimeDataHandle
  filename::String
end

function getTimeDataHandle(f::MDFFile)
  return MDFTimeDataHandle(f.filename)
end

function getindex(raw::MDFTimeDataHandle, x, y, z)
  data = h5read(raw.filename, "/measurement/dataTD", ( x, y, z) )
  return data
end
=#
