using HDF5

import HDF5: h5read

export IMTFile, IMTFileCalib, IMTFileMeas, addTrailingSingleton

abstract type IMTFile <: MPIFile end

# We use a dedicated type for calib and meas. If both versions
# are the same we use the abstract type IMTFile
mutable struct IMTFileCalib <: IMTFile
  filename::String
  file::HDF5File
end

#IMTFileCalib(filename::String, file=h5open(filename,"r")) =
#   IMTFileCalib(filename, file, nothing)

function IMTFileCalib(filename::String, file=h5open(filename,"r"))
  f = IMTFileCalib(filename, file, nothing)

 return f
end

mutable struct IMTFileMeas <: IMTFile
  filename::String
  file::HDF5File
end

function IMTFileMeas(filename::String, file=h5open(filename,"r"))
  f = IMTFileMeas(filename, file, nothing)

  return f
end

# Dispatch on file extension
function (::Type{IMTFile})(filename::String, file = h5open(filename,"r"))
 if !exists(file, "/measurements")
   return IMTFileCalib(filename, file)
 else
   return IMTFileMeas(filename, file)
 end
end


function Base.show(io::IO, f::IMTFileCalib)
  print(io, "IMT calib: ", f.filename)
end

function Base.show(io::IO, f::IMTFileMeas)
  print(io, "IMT meas: ", f.filename)
end


function getindex(f::IMTFile, parameter)
  if exists(f.file, parameter)
    return read(f.file, parameter)
  else
    return nothing
  end
end

function getindex(f::IMTFile, parameter, default)
  if exists(f.file, parameter)
    return read(f.file, parameter)
  else
    return default
  end
end


# general parameters
version(f::IMTFile) = v"0.0.0"
uuid(f::IMTFile) = uuid4()
time(f::IMTFile) = Dates.unix2datetime(0)

# study parameters
studyName(f::IMTFile) = "n.a."
studyNumber(f::IMTFile) = 0
studyUuid(f::IMTFile) = uuid4()
studyDescription(f::IMTFile) = "n.a."

# experiment parameters
experimentName(f::IMTFile) = "n.a."
experimentNumber(f::IMTFile) = 0
experimentUuid(f::IMTFile) = uuid4()
experimentDescription(f::IMTFile) = "n.a."
experimentSubject(f::IMTFile) = "n.a."
experimentIsSimulation(f::IMTFile) = true
experimentIsCalibration(f::IMTFileMeas) = false
experimentIsCalibration(f::IMTFileCalib) = true
experimentHasReconstruction(f::IMTFile) = false
experimentHasMeasurement(f::IMTFile) = true

# tracer parameters
tracerName(f::IMTFile)::Vector{String} = ["n.a"]
tracerBatch(f::IMTFile)::Vector{String} = ["n.a"]
tracerVolume(f::IMTFile)::Vector{Float64} = [0.0]
tracerConcentration(f::IMTFile)::Vector{Float64} = [0.0]
tracerSolute(f::IMTFile)::Vector{String} = ["Fe"]
#function tracerInjectionTime(f::IMTFile)::Vector{DateTime}
#  p = typeof(f) == IMTFile ? "/tracer/time" : "/tracer/injectionTime"
#  if f[p] == nothing
#    return nothing
#  end

#  if typeof(f[p]) == String
#    return [DateTime(f[p])]
#  else
#    return [DateTime(y) for y in f[p]]
#  end
#end
tracerInjectionTime(f::IMTFile) = Dates.unix2datetime(0) #DateTime( f["/tracer/injectionTime"] )
tracerVendor(f::IMTFile)::Vector{String} = ["n.a."] #[f["/tracer/vendor"]]

# scanner parameters
scannerFacility(f::IMTFile)::String = f["n.a."]
scannerOperator(f::IMTFile)::String = f["n.a."]
scannerManufacturer(f::IMTFile)::String = f["n.a."]
scannerName(f::IMTFile)::String = f["n.a."]
scannerTopology(f::IMTFile)::String = f["n.a."]

# acquisition parameters
acqStartTime(f::IMTFile)::DateTime = Dates.unix2datetime(0) #DateTime( f["/acquisition/time"] )
acqNumAverages(f::IMTFileCalib)::Int = 1 # f["/acquisition/drivefield/averages"]
acqNumAverages(f::IMTFileMeas)::Int = 1
#function acqNumFrames(f::IMTFileCalib)::Int
#  if experimentIsCalibration(f)
#    if f.mmap_measData == nothing
#      h5open(f.filename,"r") do file
#        f.mmap_measData = readmmap(file["/calibration/dataFD"])
#      end
#    end
#    return size(f.mmap_measData,2)
#  else
#    return f["/acquisition/numFrames"]
#  end
#end
acqNumFrames(f::IMTFileCalib)::Int = 1 #f["/acquisition/numFrames"]
acqNumFrames(f::IMTFileMeas)::Int = 1 #f["/acquisition/numFrames"]
acqNumPeriodsPerFrame(f::IMTFile)::Int = 1

acqGradient(f::IMTFile)::Array{Float64,4} = reshape(diagm([0.0,0.0,0.0]), 3,3,1,1)
acqOffsetField(f::IMTFile)::Array{Float64,3} = reshape([0.0,0.0,0.0],3,1,1)

# drive-field parameters
dfNumChannels(f::IMTFile) = 1
dfStrength(f::IMTFile) = 0.0 # addTrailingSingleton( addLeadingSingleton(f["/acquisition/drivefield/strength"], 2), 3)
dfPhase(f::IMTFile) = 0.0  #dfStrength(f) .*0 .+  1.5707963267948966 # Bruker specific!
dfBaseFrequency(f::IMTFile) = 2.5e6
dfCustomWaveform(f::IMTFile) = "n.a."
dfDivider(f::IMTFile) = reshape([102; 96; 99],:,1) #addTrailingSingleton(f["/acquisition/drivefield/divider"],2)
dfWaveform(f::IMTFile) = "sine"
dfCycle(f::IMTFile) = f["/timeLength"]

# receiver parameters
rxNumChannels(f::IMTFileMeas) = size(f["/measurements"],2)
rxNumChannels(f::IMTFileCalib) = 3 #size(f["/numberOfAvailableFrequencies"],2) TODO
rxBandwidth(f::IMTFile) = 1.25e6
rxNumSamplingPoints(f::IMTFile) = (f["/numberOfAvailableFrequencies"][1]-1)*2
rxTransferFunction(f::IMTFile) = nothing
rxInductionFactor(f::IMTFile) = nothing

rxUnit(f::IMTFile) = "a.u."
#rxDataConversionFactor(f::IMTFileCalib) = nothing #repeat([1.0, 0.0], outer=(1,rxNumChannels(f)))
#rxDataConversionFactor(f::IMTFileMeas) = f["/acquisition/receiver/dataConversionFactor"]

# measurements
function measData(f::IMTFileMeas, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))
  if !exists(f.file, "/measurements")
    # file is calibration
    data = f["/systemResponseFrequencies"]
    #data = f["/calibration/dataFD"]
    if ndims(data) == 4
      return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),1))
    else
      return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),size(data,5)))
    end
  end
  tdExists = exists(f.file, "/measurements")

  if tdExists
    if f.mmap_measData == nothing
      f.mmap_measData = readmmap(f.file["/measurements"])
    end
    data = zeros(Float64, rxNumSamplingPoints(f), length(receivers), length(frames))
    for (i,fr) in enumerate(frames)
      data[:,:,:,i] = f.mmap_measData[:, receivers, fr]
    end
    return reshape(data,size(data,1),size(data,2),1,size(data,3))
  else
    if f.mmap_measData == nothing
      f.mmap_measData = readmmap(f.file["/measurement/dataFD"])
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

function measData(f::IMTFileCalib, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))

  if measIsTransposed(f)
    data = f.mmap_measData[frames, :, receivers, periods]
    data = reshape(data, length(frames), size(data,2), length(receivers), length(periods))
  else
    data = f.mmap_measData[:, receivers, periods, frames]
    data = reshape(data, size(data,1), length(receivers), length(periods), length(frames))
  end
  return data
end


function measDataTDPeriods(f::IMTFileCalib, periods=1:acqNumPeriods(f),
                  receivers=1:rxNumChannels(f))
  tdExists = exists(f.file, "/measurement/dataTD")

  if tdExists
    if f.mmap_measData == nothing
      f.mmap_measData = readmmap(f.file["/measurement/dataTD"])
    end
    data = f.mmap_measData[:, receivers, periods]
    return data
  else
    if f.mmap_measData == nothing
      f.mmap_measData = readmmap(f.file["/measurement/dataFD"])
    end
    data = f.mmap_measData[:, :, receivers, periods]

    dataFD = reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4)))
    dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
    return dataTD
  end
end


function measDataTDPeriods(f::IMTFileMeas, periods=1:acqNumPeriods(f),
                  receivers=1:rxNumChannels(f))
  if measIsTransposed(f)
    error("measDataTDPeriods can currently not handle transposed data!")
  end

  data = reshape(f.mmap_measData,Val{3})[:, receivers, periods]

  return data
end

function systemMatrix(f::IMTFileCalib, rows, bgCorrection=true)
  if !experimentIsCalibration(f)
    return nothing
  end
  if f.mmap_measData == nothing
    f.mmap_measData = readmmap(f.file["/calibration/dataFD"])
  end

  data = reshape(f.mmap_measData,Val{3})[:, :, rows]
  return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3)))
end

function systemMatrix(f::IMTFileMeas, rows, bgCorrection=true)
  if !exists(f.file, "/measurement") || !measIsTransposed(f) ||
    !measIsFourierTransformed(f)
    return nothing
  end

  data_ = reshape(f.mmap_measData,size(f.mmap_measData,1),
                                  size(f.mmap_measData,2)*size(f.mmap_measData,3),
                                  size(f.mmap_measData,4))[:, rows, :]
  data = reshape(data_, Val{2})

  fgdata = data[measFGFrameIdx(f),:]
  if bgCorrection # this assumes equidistent bg frames
    println("Applying bg correction on system matrix (MDF)")
    bgdata = data[measBGFrameIdx(f),:]
    bgdataInterp = interpolate(bgdata, (BSpline(Linear()), NoInterp()))
    #Cubic does not work for complex numbers
    origIndex = measFramePermutation(f)
    M = size(fgdata,1)
    K = size(bgdata,1)
    N = M + K
    for m=1:M
      alpha = (origIndex[m]-1)/(N-1)*(K-1)+1
      for k=1:size(fgdata,2)
        fgdata[m,k] -= bgdataInterp(alpha,k)
      end
    end
  end
  return fgdata
end

function systemMatrixWithBG(f::IMTFileMeas)
  if !exists(f.file, "/measurement") || !measIsTransposed(f) ||
      !measIsFourierTransformed(f)
      return nothing
  end

  data = f.mmap_measData[:, :, :, :]
  return data
end


function measIsFourierTransformed(f::IMTFileCalib)
  if !experimentIsCalibration(f)
    return false
  else
    return true
  end
end

measIsFourierTransformed(f::IMTFile) = true
measIsTFCorrected(f::IMTFile) = false
measIsSpectralLeakageCorrected(f::IMTFile) = false

measIsBGCorrected(f::IMTFileCalib) = false
measIsBGCorrected(f::IMTFileMeas) = true #Bool(f["/measurement/isBackgroundCorrected"])

measIsFrequencySelection(f::IMTFile) = false

measIsTransposed(f::IMTFileCalib) = true
measIsTransposed(f::IMTFileMeas) = false

measIsFramePermutation(f::IMTFileCalib) = true
measIsFramePermutation(f::IMTFileMeas) = false

measIsBGFrame(f::IMTFile) = zeros(Bool, acqNumFrames(f))

measFramePermutation(f::IMTFileCalib) = nothing
measFramePermutation(f::IMTFileMeas) = nothing

#fullFramePermutation(f::IMTFile) = fullFramePermutation(f, calibIsMeanderingGrid(f))

#calibrations
#calibSNR(f::IMTFile) = addTrailingSingleton(f["/calibration/snrFD"],3)
calibFov(f::IMTFile) = f["/fov"]
calibFovCenter(f::IMTFile) = [0.0,0.0,0.0]
calibSize(f::IMTFile) = nothing
calibOrder(f::IMTFile) = "xyz"
#calibOffsetField(f::IMTFile) = [0.0,0.0,0.0]
calibDeltaSampleSize(f::IMTFile) = [0.0,0.0,0.0]
calibMethod(f::IMTFile) = "simulation"
#calibIsMeanderingGrid(f::IMTFile) = Bool(f["/calibration/isMeanderingGrid", 0])
#calibPositions(f::IMTFile) = f["/calibration/positions"]


# additional functions that should be implemented by an MPIFile
filepath(f::IMTFile) = f.filename
