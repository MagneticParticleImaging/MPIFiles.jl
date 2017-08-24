using HDF5

export MDFFile, MDFFileV1, MDFFileV2, addTrailingSingleton, addLeadingSingleton

@compat abstract type MDFFile <: MPIFile end

# We use a dedicated type for v1 and v2. If both versions
# are the same we use the abstract type MDFFile
type MDFFileV1 <: MDFFile
  filename::String
  param_cache
  mmap_measData
end

MDFFileV1(filename::String) = MDFFileV1(filename,Dict{String,Any}(),nothing)

type MDFFileV2 <: MDFFile
  filename::String
  param_cache
  mmap_measData
end

MDFFileV2(filename::String) = MDFFileV2(filename,Dict{String,Any}(),nothing)

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
uuid(f::MDFFile) = str2uuid(f["/uuid"])
time(f::MDFFileV1) = DateTime( f["/date"] )
time(f::MDFFileV2) = DateTime( f["/time"] )

# study parameters
studyName(f::MDFFile) = f["/study/name"]
studyNumber(f::MDFFileV1) = 0
studyNumber(f::MDFFileV2) = f["/study/number"]
studyUuid(f::MDFFileV1) = nothing
studyUuid(f::MDFFileV2) = str2uuid(f["/study/uuid"])
studyDescription(f::MDFFileV1) = "n.a."
studyDescription(f::MDFFileV2) = f["/study/description"]

# experiment parameters
experimentName(f::MDFFileV1) = "n.a."
experimentName(f::MDFFileV2) = f["/experiment/name"]
experimentNumber(f::MDFFileV1) = parse(Int64, f["/study/experiment"])
experimentNumber(f::MDFFileV2) = f["/experiment/number"]
experimentUuid(f::MDFFileV1) = nothing
experimentUuid(f::MDFFileV2) = str2uuid(f["/experiment/uuid"])
experimentDescription(f::MDFFileV1) = f["/study/description"]
experimentDescription(f::MDFFileV2) = f["/experiment/description"]
experimentSubject(f::MDFFileV1) = f["/study/subject"]
experimentSubject(f::MDFFileV2) = f["/experiment/subject"]
experimentIsSimulation(f::MDFFileV2) = Bool( f["/experiment/isSimulation"] )
experimentIsSimulation(f::MDFFileV1) = Bool( f["/study/simulation"] )
experimentIsCalibration(f::MDFFile) = h5exists(f.filename, "/calibration")
experimentHasReconstruction(f::MDFFile) = h5exists(f.filename, "/reconstruction")
experimentHasMeasurement(f::MDFFileV1) = h5exists(f.filename, "/measurement") ||
                                         h5exists(f.filename, "/calibration")
experimentHasMeasurement(f::MDFFileV2) = h5exists(f.filename, "/measurement")

_makeStringArray(s::String) = [s]
_makeStringArray{T<:AbstractString}(s::Vector{T}) = s

# tracer parameters
tracerName(f::MDFFileV1) = [f["/tracer/name"]]
tracerName(f::MDFFileV2) = _makeStringArray(f["/tracer/name"])
tracerBatch(f::MDFFileV1) = [f["/tracer/batch"]]
tracerBatch(f::MDFFileV2) = _makeStringArray(f["/tracer/batch"])
tracerVolume(f::MDFFileV1) = [f["/tracer/volume"]]
tracerVolume(f::MDFFileV2) = [f["/tracer/volume"]...]
tracerConcentration(f::MDFFileV1) = [f["/tracer/concentration"]]
tracerConcentration(f::MDFFileV2) = [f["/tracer/concentration"]...]
tracerSolute(f::MDFFileV2) = _makeStringArray(f["/tracer/solute"])
tracerSolute(f::MDFFileV1) = ["Fe"]
function tracerInjectionTime(f::MDFFile)
  p = typeof(f) == MDFFileV1 ? "/tracer/time" : "/tracer/injectionTime"
  if f[p] == nothing
    return nothing
  end

  if typeof(f[p]) == String
    return [DateTime(f[p])]
  else
    return [DateTime(y) for y in f[p]]
  end
end
#tracerInjectionTime(f::MDFFileV2) = DateTime( f["/tracer/injectionTime"] )
tracerVendor(f::MDFFileV1) = [f["/tracer/vendor"]]
tracerVendor(f::MDFFileV2) = _makeStringArray(f["/tracer/vendor"])

# scanner parameters
scannerFacility(f::MDFFile) = f["/scanner/facility"]
scannerOperator(f::MDFFile) = f["/scanner/operator"]
scannerManufacturer(f::MDFFile) = f["/scanner/manufacturer"]
scannerName(f::MDFFileV1) = f["/scanner/model"]
scannerName(f::MDFFileV2) = f["/scanner/name"]
scannerTopology(f::MDFFile) = f["/scanner/topology"]

# acquisition parameters
acqStartTime(f::MDFFileV1) = DateTime( f["/acquisition/time"] )
acqStartTime(f::MDFFileV2) = DateTime( f["/acquisition/startTime"] )
acqFramePeriod(f::MDFFile) = f["/acquisition/framePeriod"]
acqNumPatches(f::MDFFile) = f["/acquisition/numPatches"]
acqNumAverages(f::MDFFileV1) = f["/acquisition/drivefield/averages"]
acqNumAverages(f::MDFFileV2) = f["/acquisition/numAverages"]
function acqNumFrames(f::MDFFileV1)
  if experimentIsCalibration(f)
    if f.mmap_measData == nothing
      h5open(f.filename,"r") do file
        f.mmap_measData = readmmap(file["/calibration/dataFD"])
      end
    end
    return size(f.mmap_measData,2)
  else
    return f["/acquisition/numFrames"]
  end
end
acqNumFrames(f::MDFFileV2) = f["/acquisition/numFrames"]
acqNumPeriods(f::MDFFileV1) = 1
acqNumPeriods(f::MDFFileV2) = f["/acquisition/numPeriods"]

acqGradient(f::MDFFileV1) = addTrailingSingleton(f["/acquisition/gradient"],2)
acqGradient(f::MDFFileV2) = f["/acquisition/gradient"]
acqOffsetField(f::MDFFile) = f["/acquisition/offsetField"]
acqOffsetFieldShift(f::MDFFileV1) = addTrailingSingleton(
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
rxBandwidth(f::MDFFile) = f["/acquisition/receiver/bandwidth"]
rxNumSamplingPoints(f::MDFFile) = f["/acquisition/receiver/numSamplingPoints"]
function rxTransferFunction(f::MDFFile)
  parameter = "/acquisition/receiver/transferFunction"
  if h5exists(f.filename, parameter)
    return readComplexArray(f.filename, parameter)
  else
    return nothing
  end
end
rxInductionFactor(f::MDFFileV1) = nothing
rxInductionFactor(f::MDFFileV2) = f["/acquisition/receiver/inductionFactor"]

rxUnit(f::MDFFileV1) = "a.u."
rxUnit(f::MDFFileV2) = f["/acquisition/receiver/unit"]
rxDataConversionFactor(f::MDFFileV1) = repeat([1.0, 0.0], outer=(1,rxNumChannels(f)))
rxDataConversionFactor(f::MDFFileV2) = f["/acquisition/receiver/dataConversionFactor"]

# measurements
function measData(f::MDFFileV1, frames=1:acqNumFrames(f), periods=1:acqNumPeriods(f),
                  receivers=1:rxNumChannels(f))
  if !h5exists(f.filename, "/measurement")
    # the V1 file is a calibration
    data = f["/calibration/dataFD"]
    if ndims(data) == 4
      return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),1))
    else
      return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),size(data,5)))
    end
  end
  tdExists = h5exists(f.filename, "/measurement/dataTD")

  if tdExists
    if f.mmap_measData == nothing
      h5open(f.filename,"r") do file
        f.mmap_measData = readmmap(file["/measurement/dataTD"])
      end
    end
    data = zeros(Float64, rxNumSamplingPoints(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:, receivers, fr]
    end
    return reshape(data,size(data,1),size(data,2),1,size(data,3))
  else
    if f.mmap_measData == nothing
      h5open(f.filename,"r") do file
        f.mmap_measData = readmmap(file["/measurement/dataFD"])
      end
    end
    data = zeros(Float64, 2, rxNumFrequencies(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:,:,receivers, fr]
    end

    dataFD = reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4)))
    dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
    return reshape(dataTD,size(dataTD,1),size(dataTD,2),1,size(dataTD,3))
  end
end

function measData(f::MDFFileV2, frames=1:acqNumFrames(f), periods=1:acqNumPeriods(f),
                  receivers=1:rxNumChannels(f))
  if !h5exists(f.filename, "/measurement")
    return nothing
  end
  if f.mmap_measData == nothing
    h5open(f.filename,"r") do file
      parameter = "/measurement/data"
      if !isComplexArray(file, parameter)
        f.mmap_measData = readmmap(file[parameter])
      else
        f.mmap_measData = readmmap(file[parameter], Array{getComplexType(file,parameter)} )
      end
    end
  end

  if measIsTransposed(f)
    data = f.mmap_measData[frames, :, receivers, patches]
    data = reshape(data, length(frames), size(data,2), length(receivers), length(periods))
  else
    data = f.mmap_measData[:, receivers, patches, frames]
    data = reshape(data, size(data,1), length(receivers), length(periods), length(frames))
  end
  return data
end

function systemMatrix(f::MDFFileV1, rows, bgCorrection=true)
  if !experimentIsCalibration(f)
    return nothing
  end
  if f.mmap_measData == nothing
    h5open(f.filename,"r") do file
      f.mmap_measData = readmmap(file["/calibration/dataFD"])
    end
  end

  data = reshape(f.mmap_measData,Val{3})[:, :, rows]
  return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3)))
end

function systemMatrix(f::MDFFileV2, rows, bgCorrection=true)
  if !h5exists(f.filename, "/measurement") || !measIsTransposed(f) ||
    !measIsFourierTransformed(f)
    return nothing
  end
  if f.mmap_measData == nothing
    h5open(f.filename,"r") do file
      parameter = "/measurement/data"
      if !isComplexArray(file, parameter)
        f.mmap_measData = readmmap(file[parameter])
      else
        f.mmap_measData = readmmap(file[parameter], Array{getComplexType(file,parameter)} )
      end
    end
  end
  data = reshape(f.mmap_measData,Val{2})[:, rows]

  fgdata = data[measFGFrameIdx(f),:]
  if bgCorrection
    bgdata = data[measBGFrameIdx(f),:]
    fgdata[:,:] .-= mean(bgdata,1)
  end
  return fgdata
end

function systemMatrixWithBG(f::MDFFileV2)
  if !h5exists(f.filename, "/measurement") || !measIsTransposed(f) ||
      !measIsFourierTransformed(f)
      return nothing
  end
  if f.mmap_measData == nothing
    h5open(f.filename,"r") do file
      parameter = "/measurement/data"
      if !isComplexArray(file, parameter)
        f.mmap_measData = readmmap(file[parameter])
      else
        f.mmap_measData = readmmap(file[parameter], Array{getComplexType(file,parameter)} )
      end
    end
  end

  data = f.mmap_measData[:, :, :, :]
  return data
end


function measIsFourierTransformed(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsFourierTransformed(f::MDFFileV2) = Bool(f["/measurement/isFourierTransformed"])

measIsTFCorrected(f::MDFFileV1) = false
measIsTFCorrected(f::MDFFileV2) = Bool(f["/measurement/isTransferFunctionCorrected"])

measIsSpectralLeakageCorrected(f::MDFFileV1) = false
measIsSpectralLeakageCorrected(f::MDFFileV2) = Bool(f["/measurement/isSpectralLeakageCorrected"])

function measIsBGCorrected(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsBGCorrected(f::MDFFileV2) = Bool(f["/measurement/isBackgroundCorrected"])

measIsFrequencySelection(f::MDFFileV1) = false
measIsFrequencySelection(f::MDFFileV2) = Bool(f["/measurement/isFrequencySelection"])

function measIsTransposed(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsTransposed(f::MDFFileV2) = Bool(f["/measurement/isTransposed"])

function measIsFramePermutation(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsFramePermutation(f::MDFFileV2) = f["/measurement/isFramePermutation"]
measIsBGFrame(f::MDFFileV1) = zeros(Bool, acqNumFrames(f))
measIsBGFrame(f::MDFFileV2) = convert(Array{Bool},f["/measurement/isBackgroundFrame"])
measFramePermutation(f::MDFFileV1) = nothing
measFramePermutation(f::MDFFileV2) = f["/measurement/framePermutation"]

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


