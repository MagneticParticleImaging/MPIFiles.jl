export MDFFile, MDFFileV1, MDFFileV2, addTrailingSingleton, addLeadingSingleton

abstract type MDFFile <: MPIFile end

# We use a dedicated type for v1 and v2. If both versions
# are the same we use the abstract type MDFFile
mutable struct MDFFileV1{T} <: MDFFile
  filename::String
  file::HDF5.File
  mmap_measData::T
end

function MDFFileV1(filename::String, file=h5open(filename,"r"))
  tdExists = haskey(file, "/measurement/dataTD")
  if haskey(file, "/measurement/dataTD")
    mmap_measData = HDF5.readmmap(file["/measurement/dataTD"])
  elseif haskey(file, "/measurement/dataFD")
    mmap_measData = HDF5.readmmap(file["/measurement/dataFD"])
  elseif haskey(file, "/calibration/dataFD")
    mmap_measData = HDF5.readmmap(file["/calibration/dataFD"])
  else
    mmap_measData = nothing
  end

  f = MDFFileV1(filename, file, mmap_measData)
end

mutable struct MDFFileV2{T} <: MDFFile
  filename::String
  file::HDF5.File
  mmap_measData::T
end

function MDFFileV2(filename::String, file=h5open(filename,"r"))
  parameter = "/measurement/data"
  if haskey(file, "/measurement/data")
    if !isComplexArray(file, parameter)
      mmap_measData = HDF5.readmmap(file[parameter])
    else
      #mmap_measData = HDF5.readmmap(file[parameter], Array{getComplexType(file,parameter)} )
      mmap_measData = HDF5.readmmap(file[parameter], getComplexType(file,parameter) )
    end
  else
    mmap_measData = nothing
  end

  f = MDFFileV2(filename, file, mmap_measData)

  return f
end

# This dispatches on the file extension and automatically
# generates the correct type
function MDFFile(filename::String, file = h5open(filename,"r"))
  vers = VersionNumber( read(file, "/version") )
  if vers < v"2.0"
    return MDFFileV1(filename, file)
  else
    return MDFFileV2(filename, file)
  end
end

Base.close(f::MDFFile) = close(f.file)


function h5haskey(filename, parameter)
  return h5open(filename) do file
    haskey(file, parameter)
  end
end

function getindex(f::MDFFile, parameter)
  if haskey(f.file, parameter)
    return read(f.file, parameter)
  else
    return nothing
  end
end

function getindex(f::MDFFile, parameter, default)
  #if !haskey(f.param_cache,parameter)
  #  f.param_cache[parameter] = h5read_(f.filename, parameter, default)
  #end
  #return f.param_cache[parameter]
  if haskey(f.file, parameter)
    return read(f.file, parameter)
  else
    return default
  end
end



# general parameters
version(f::MDFFile)::Union{VersionNumber, Missing} = @keyrequired VersionNumber( f["/version"] )
uuid(f::MDFFile)::Union{UUID, Missing} = @keyrequired UUID(f["/uuid"])
time(f::MDFFileV1)::Union{DateTime, Missing} = @keyrequired DateTime( f["/date"] )
time(f::MDFFileV2)::Union{DateTime, Missing} = @keyrequired DateTime( f["/time"] )

# study parameters
studyName(f::MDFFile)::Union{String, Missing} = @keyrequired f["/study/name"]
studyNumber(f::MDFFileV1)::Union{Int, Missing} = 0
studyNumber(f::MDFFileV2)::Union{Int, Missing} = @keyrequired f["/study/number"]
studyUuid(f::MDFFileV1) = nothing
studyUuid(f::MDFFileV2) = @keyrequired UUID(f["/study/uuid"])
studyDescription(f::MDFFileV1)::Union{String, Missing} = "n.a."
studyDescription(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/study/description"]
function studyTime(f::MDFFile)
  t = f["/study/time"]
  if typeof(t)==String
   return DateTime(t)
  else
   return nothing
  end
end

# experiment parameters
experimentName(f::MDFFileV1)::Union{String, Missing} = "n.a."
experimentName(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/experiment/name"]
experimentNumber(f::MDFFileV1)::Union{Int64, Missing} = @keyrequired parse(Int64, f["/study/experiment"])
experimentNumber(f::MDFFileV2)::Union{Int64, Missing} = @keyrequired f["/experiment/number"]
experimentUuid(f::MDFFileV1) = nothing
experimentUuid(f::MDFFileV2) = @keyrequired UUID(f["/experiment/uuid"])
experimentDescription(f::MDFFileV1)::Union{String, Missing} = @keyrequired f["/study/description"]
experimentDescription(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/experiment/description"]
experimentSubject(f::MDFFileV1)::Union{String, Missing} = @keyrequired f["/study/subject"]
experimentSubject(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/experiment/subject"]
experimentIsSimulation(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool( f["/experiment/isSimulation"] )
experimentIsSimulation(f::MDFFileV1)::Union{Bool, Missing} = @keyrequired Bool( f["/study/simulation"] )
experimentIsCalibration(f::MDFFile)::Bool = haskey(f.file, "/calibration") && haskey(f.file, "/calibration/method") # Last part is a little workaround for the testcases, since the file somehow contains a `deltaSampleSize` field
experimentHasReconstruction(f::MDFFile)::Bool = haskey(f.file, "/reconstruction")
experimentHasMeasurement(f::MDFFileV1)::Bool = haskey(f.file, "/measurement") ||
                                         haskey(f.file, "/calibration")
experimentHasMeasurement(f::MDFFileV2)::Bool = haskey(f.file, "/measurement")

_makeStringArray(s::String) = [s]
_makeStringArray(s::Vector{T}) where {T<:AbstractString} = s

# tracer parameters
tracerName(f::MDFFileV1)::Union{Vector{String}, Missing} = @keyrequired [f["/tracer/name"]]
tracerName(f::MDFFileV2)::Union{Vector{String}, Missing} = @keyrequired _makeStringArray(f["/tracer/name"])
tracerBatch(f::MDFFileV1)::Union{Vector{String}, Missing} = @keyrequired [f["/tracer/batch"]]
tracerBatch(f::MDFFileV2)::Union{Vector{String}, Missing} = @keyrequired _makeStringArray(f["/tracer/batch"])
tracerVolume(f::MDFFileV1)::Union{Vector{Float64}, Missing} = @keyrequired [f["/tracer/volume"]]
tracerVolume(f::MDFFileV2)::Union{Vector{Float64}, Missing} = @keyrequired [f["/tracer/volume"]...]
tracerConcentration(f::MDFFileV1)::Union{Vector{Float64}, Missing} = @keyrequired [f["/tracer/concentration"]]
tracerConcentration(f::MDFFileV2)::Union{Vector{Float64}, Missing} = @keyrequired [f["/tracer/concentration"]...]
tracerSolute(f::MDFFileV2)::Union{Vector{String}, Missing} = @keyrequired _makeStringArray(f["/tracer/solute"])
tracerSolute(f::MDFFileV1)::Union{Vector{String}, Missing} = ["Fe"]
function tracerInjectionTime(f::MDFFile)::Union{Vector{DateTime}, Nothing}
  p = typeof(f) <: MDFFileV1 ? "/tracer/time" : "/tracer/injectionTime"
  if isnothing(f[p])
    return nothing
  end

  if typeof(f[p]) == String
    return [DateTime(f[p])]
  else
    return [DateTime(y) for y in f[p]]
  end
end
#tracerInjectionTime(f::MDFFileV2) = DateTime( f["/tracer/injectionTime"] )
tracerVendor(f::MDFFileV1)::Union{Vector{String}, Missing} = @keyrequired [f["/tracer/vendor"]]
tracerVendor(f::MDFFileV2)::Union{Vector{String}, Missing} = @keyrequired _makeStringArray(f["/tracer/vendor"])

# scanner parameters
scannerBoreSize(f::MDFFile)::Union{String, Nothing} = @keyoptional f["/scanner/boreSize"]
scannerFacility(f::MDFFile)::Union{String, Missing} = @keyrequired f["/scanner/facility"]
scannerOperator(f::MDFFile)::Union{String, Missing} = @keyrequired f["/scanner/operator"]
scannerManufacturer(f::MDFFile)::Union{String, Missing} = @keyrequired f["/scanner/manufacturer"]
scannerName(f::MDFFileV1)::Union{String, Missing} = @keyrequired f["/scanner/model"]
scannerName(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/scanner/name", ""]
scannerTopology(f::MDFFile)::Union{String, Missing} = @keyrequired f["/scanner/topology"]

# acquisition parameters
acqStartTime(f::MDFFileV1)::Union{DateTime, Nothing} = @keyoptional DateTime( f["/acquisition/time"] )
acqStartTime(f::MDFFileV2)::Union{DateTime, Nothing} = @keyoptional DateTime( f["/acquisition/startTime"] )
acqNumAverages(f::MDFFileV1)::Union{Int, Missing} = @keyrequired f["/acquisition/drivefield/averages"]
acqNumAverages(f::MDFFileV2)::Union{Int, Missing} = @keyrequired f["/acquisition/numAverages",1]
function acqNumFrames(f::MDFFileV1)::Int
  if experimentIsCalibration(f)
    return size(f.mmap_measData,2)
  else
    return f["/acquisition/numFrames"]
  end
end
acqNumFrames(f::MDFFileV2)::Union{Int, Missing} = @keyrequired f["/acquisition/numFrames"]
acqNumPeriodsPerFrame(f::MDFFileV1)::Union{Int, Missing} = 1
function acqNumPeriodsPerFrame(f::MDFFileV2)::Union{Int, Missing}
  if haskey(f.file, "/acquisition/numPeriodsPerFrame")
    return f["/acquisition/numPeriodsPerFrame"]
  elseif haskey(f.file, "/acquisition/numPeriods")
    return f["/acquisition/numPeriods"]
  else
    return 1
  end
end

acqGradient(f::MDFFileV1)::Union{Array{Float64,4}, Nothing} = @keyoptional reshape(Matrix(Diagonal(f["/acquisition/gradient"])), 3,3,1,1)
function acqGradient(f::MDFFileV2)::Union{Array{Float64,4}, Nothing}
  G = f["/acquisition/gradient", Matrix(Diagonal([1.0,1.0,-2.0]))]
  if ndims(G) == 4
   return G
  elseif ndims(G) == 3 # for corrupt files
   return reshape(G,3,3,1,size(G,3))
  elseif ndims(G) == 2 && prod(size(G)) == 9  # for corrupt files
   return reshape(G,3,3,1,1)
  else # for corrupt files
   return reshape(Matrix(Diagonal(vec(G))),3,3,1,1)
  end
end

acqOffsetField(f::MDFFileV1)::Union{Array{Float64,3}, Nothing} = @keyoptional f["/acquisition/offsetField", reshape([0.0,0.0,0.0],3,1,1)  ]
function acqOffsetField(f::MDFFileV2)::Union{Array{Float64,3}, Nothing}
  H = f["/acquisition/offsetField", reshape([0.0,0.0,0.0],3,1,1)  ]
  if ndims(H) == 3
   return H
  else # for corrupt files
   return reshape(H,:,1,1)
  end
end

# drive-field parameters
dfNumChannels(f::MDFFile)::Union{Int, Missing} = @keyrequired f["/acquisition/drivefield/numChannels"]
dfStrength(f::MDFFileV1)::Union{Array{Float64,3}, Missing} = @keyrequired addTrailingSingleton( addLeadingSingleton(
         f["/acquisition/drivefield/strength"], 2), 3)
dfStrength(f::MDFFileV2)::Union{Array{Float64,3}, Missing} = @keyrequired f["/acquisition/drivefield/strength"]
dfPhase(f::MDFFileV1)::Union{Array{Float64,3}, Missing} = dfStrength(f) .*0 .+  1.5707963267948966 # Bruker specific!
dfPhase(f::MDFFileV2)::Union{Array{Float64,3}, Missing} = @keyrequired f["/acquisition/drivefield/phase"]
dfBaseFrequency(f::MDFFile)::Union{Float64, Missing} = @keyrequired f["/acquisition/drivefield/baseFrequency"]
dfCustomWaveform(f::MDFFileV2)::Union{String, Nothing} = @keyoptional f["/acquisition/drivefield/customWaveform"]
dfDivider(f::MDFFileV1) = @keyrequired addTrailingSingleton(
                f["/acquisition/drivefield/divider"],2)
dfDivider(f::MDFFileV2) = f["/acquisition/drivefield/divider"]
dfWaveform(f::MDFFileV1) = "sine"
function dfWaveform(f::MDFFileV2)::Union{Array{String, 2}, Missing}
  value = @keyrequired f["/acquisition/drivefield/waveform"]
  if typeof(value) == String # Fix wrong implementation of specs downstream!
    return fill(value, (1, 1))
  else
    return value # This should be an 2D array of strings or missing
  end
end

dfCycle(f::MDFFile)::Union{Float64, Missing} = @keyrequired f["/acquisition/drivefield/cycle"]
dfCycle(f::MDFFileV1)::Union{Float64, Missing} = @keyrequired f["/acquisition/drivefield/period"]

# receiver parameters
rxNumChannels(f::MDFFile)::Union{Int64, Missing} = @keyrequired f["/acquisition/receiver/numChannels"]
rxBandwidth(f::MDFFile)::Union{Float64, Missing} = @keyrequired f["/acquisition/receiver/bandwidth"]
rxNumSamplingPoints(f::MDFFile)::Union{Int64, Missing}= @keyrequired f["/acquisition/receiver/numSamplingPoints"]
function rxTransferFunction(f::MDFFile)
  parameter = "/acquisition/receiver/transferFunction"
  if haskey(f.file, parameter)
    return readComplexArray(f.filename, parameter)
  else
    return nothing
  end
end
function rxTransferFunctionFileName(f::MDFFile)
  parameter = "/acquisition/receiver/transferFunctionFileName"
  if haskey(f.file, parameter)
    return f[parameter]
  else
    return nothing
  end
end
function rxHasTransferFunction(f::MDFFile)
  haskey(f.file, "/acquisition/receiver/transferFunction")
end
rxInductionFactor(f::MDFFileV1) = nothing
rxInductionFactor(f::MDFFileV2) = @keyrequired f["/acquisition/receiver/inductionFactor"]

rxUnit(f::MDFFileV1)::Union{String, Missing} = "a.u."
rxUnit(f::MDFFileV2)::Union{String, Missing} = @keyrequired f["/acquisition/receiver/unit"]
rxDataConversionFactor(f::MDFFileV1) = repeat([1.0, 0.0], outer=(1,rxNumChannels(f)))
rxDataConversionFactor(f::MDFFileV2) = @keyrequired f["/acquisition/receiver/dataConversionFactor"]

# measurements
function measData(f::MDFFileV1, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))
  if !haskey(f.file, "/measurement")
    # the V1 file is a calibration
    data = f["/calibration/dataFD"]
    if ndims(data) == 4
      return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4),1))
    else
      return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4),size(data,5)))
    end
  end
  tdExists = haskey(f.file, "/measurement/dataTD")

  if tdExists
    data = zeros(Float64, rxNumSamplingPoints(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:, receivers, fr]
    end
    return reshape(data,size(data,1),size(data,2),1,size(data,3))
  else
    data = zeros(Float64, 2, rxNumFrequencies(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:,:,receivers, fr]
    end

    dataFD = reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4)))
    dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
    return reshape(dataTD,size(dataTD,1),size(dataTD,2),1,size(dataTD,3))
  end
end

measDataRaw(f::MDFFileV2) = f.mmap_measData

# Functions working on both MDFFilev2 and MDFv2InMemory were moved to MDFCommon.jl

function measDataTDPeriods(f::MDFFileV1, periods=1:acqNumPeriods(f),
  receivers=1:rxNumChannels(f))
  tdExists = haskey(f.file, "/measurement/dataTD")

  if tdExists
  data = f.mmap_measData[:, receivers, periods]
  return data
  else
  data = f.mmap_measData[:, :, receivers, periods]

  dataFD = reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3),size(data,4)))
  dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
  return dataTD
  end
end

function systemMatrix(f::MDFFileV1, rows, bgCorrection=true)
  if !experimentIsCalibration(f)
    return nothing
  end

  data = reshape(f.mmap_measData,Val(3))[:, :, rows]
  return reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data,2),size(data,3)))
end


function measIsFourierTransformed(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsFourierTransformed(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isFourierTransformed"])

measIsTFCorrected(f::MDFFileV1)::Union{Bool, Missing} = false
measIsTFCorrected(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isTransferFunctionCorrected"])

measIsSpectralLeakageCorrected(f::MDFFileV1)::Union{Bool, Missing} = @keyrequired false
measIsSpectralLeakageCorrected(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isSpectralLeakageCorrected"])

function measIsBGCorrected(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsBGCorrected(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isBackgroundCorrected"])

measIsFrequencySelection(f::MDFFileV1)::Union{Bool, Missing} = false
measIsFrequencySelection(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isFrequencySelection"])
measFrequencySelection(f::MDFFileV2)::Union{Vector{Int64}, Nothing} = @keyoptional f["/measurement/frequencySelection"]

measIsSparsityTransformed(f::MDFFileV1) = false
function measIsSparsityTransformed(f::MDFFileV2)
  if haskey(f.file, "/measurement/isSparsityTransformed")
    Bool(f["/measurement/isSparsityTransformed"])
  else
    return false
  end
end

function measIsFastFrameAxis(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end

function measIsFastFrameAxis(f::MDFFileV2)
  if haskey(f.file, "/measurement/isFastFrameAxis")
    return Bool(f["/measurement/isFastFrameAxis"])
  else
    @warn "/measurement/isFastFrameAxis missing in MDF data set. `measIsFastFrameAxis` returning false per default."
    return false
  end
end

function measIsFramePermutation(f::MDFFileV1)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end
measIsFramePermutation(f::MDFFileV2)::Union{Bool, Missing} = @keyrequired Bool(f["/measurement/isFramePermutation"])
measIsBGFrame(f::MDFFileV1)::Union{Vector{Bool}, Missing} = zeros(Bool, acqNumFrames(f))
measIsBGFrame(f::MDFFileV2)::Union{Vector{Bool}, Missing} = @keyrequired convert(Array{Bool},f["/measurement/isBackgroundFrame"])
measFramePermutation(f::MDFFileV1)::Union{Vector{Int64}, Nothing} = nothing
measFramePermutation(f::MDFFileV2)::Union{Vector{Int64}, Nothing} = @keyoptional f["/measurement/framePermutation"]
measSparsityTransformation(f::MDFFileV1)::Union{String, Nothing} = nothing
measSparsityTransformation(f::MDFFileV2)::Union{String, Nothing} = @keyoptional f["/measurement/sparsityTransformation"]
measSubsamplingIndices(f::MDFFileV2)::Union{Array{Integer, 4}, Nothing} = @keyoptional f["/measurement/subsamplingIndices"]

fullFramePermutation(f::MDFFile) = fullFramePermutation(f, calibIsMeanderingGrid(f))

measIsCalibProcessed(f::MDFFile) = measIsFramePermutation(f) && 
                                   measIsFourierTransformed(f) &&
                                   measIsFastFrameAxis(f)
measTemperatures(f::MDFFile) = @keyoptional f["/measurement/_temperatures"] # non-standard

#calibrations
calibSNR(f::MDFFileV1) = @keyoptional addTrailingSingleton(f["/calibration/snrFD"],3)
calibSNR(f::MDFFileV2) = @keyoptional f["/calibration/snr"]
calibFov(f::MDFFile) = @keyoptional f["/calibration/fieldOfView"]
calibFovCenter(f::MDFFile) = @keyoptional f["/calibration/fieldOfViewCenter"]
calibSize(f::MDFFile) = @keyoptional f["/calibration/size"]
calibOrder(f::MDFFile) = @keyoptional f["/calibration/order"]
calibOffsetField(f::MDFFile) = @keyoptional f["/calibration/offsetField"]
calibDeltaSampleSize(f::MDFFile) = @keyoptional f["/calibration/deltaSampleSize",[0.0,0.0,0.0]]
calibMethod(f::MDFFile) = @keyrequired f["/calibration/method"]
calibIsMeanderingGrid(f::MDFFile) = @keyoptional Bool(f["/calibration/isMeanderingGrid", 0])
calibPositions(f::MDFFile) = @keyoptional f["/calibration/positions"]

# reconstruction results
recoData(f::MDFFileV1) = @keyrequired addLeadingSingleton(f[ "/reconstruction/data"], 3)
recoData(f::MDFFileV2) = @keyrequired f["/reconstruction/data"]
recoFov(f::MDFFile)::Union{Vector{Float64}, Nothing} = @keyoptional f["/reconstruction/fieldOfView"]
recoFovCenter(f::MDFFile)::Union{Vector{Float64}, Nothing} = @keyoptional f["/reconstruction/fieldOfViewCenter"]
recoSize(f::MDFFile)::Union{Vector{Int64}, Nothing} = @keyoptional f["/reconstruction/size"]
recoOrder(f::MDFFile)::Union{String, Nothing} = @keyoptional f["/reconstruction/order"]
recoPositions(f::MDFFile)::Union{Array{Float64, 2}, Nothing} = @keyoptional f["/reconstruction/positions"]
recoIsOverscanRegion(f::MDFFileV1)::Union{Vector{Bool}, Nothing} = nothing
recoIsOverscanRegion(f::MDFFileV2)::Union{Vector{Bool}, Nothing} = @keyoptional f["/reconstruction/isOverscanRegion"]

# this is non-standard
function recoParameters(f::MDFFile)
  if !haskey(f.file, "/reconstruction/_parameters")
    return nothing
  else
    return loadParams(f.file, "/reconstruction/_parameters")
  end
end

auxiliaryData(f::MDFFileV2) = @keyoptional f["/custom/auxiliaryData"]

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
