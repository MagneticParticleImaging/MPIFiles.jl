using HDF5

export MDFFile, MDFFileV1, MDFFileV2

abstract MDFFile <: MPIFile

# We use a dedicated type for v1 and v2. If both versions
# are the same we use the abstract type MDFFile
type MDFFileV1 <: MDFFile
  filename::String
end

type MDFFileV2 <: MDFFile
  filename::String
end

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
  print(io, "MDF v1: ", f.path)
end

function Base.show(io::IO, f::MDFFileV2)
  print(io, "MDF v2: ", f.path)
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
   return h5readornull(f.filename, parameter)
end



# general parameters
version(f::MDFFile) = VersionNumber( f["/version"] )
uuid(f::MDFFile) = f["/uuid"]
time(f::MDFFile) = DateTime( f["/date"] )

# study parameters
studyName(f::MDFFile) = f["/study/name"]
studyExperiment(f::MDFFileV1) = parse(Int64, f["/study/experiment"])
studyExperiment(f::MDFFileV2) = f["/study/experiment"]
studyDescription(f::MDFFile) = f["/study/description"]
studySubject(f::MDFFile) = f["/study/subject"]
studyIsSimulation(f::MDFFileV2) = Bool( f["/study/isSimulation"] )
studyIsSimulation(f::MDFFileV1) = Bool( f["/study/simulation"] )
studyIsCalibration(f::MDFFileV2) = Bool( f["/study/isCalibration"] )
studyIsCalibration(f::MDFFileV1) = h5exists(f.filename, "/calibration")

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
acqFov(f::MDFFileV1) = addLeadingSingleton( f["/acquisition/drivefield/fieldOfView"],2 )
acqFov(f::MDFFileV2) = f["/acquisition/fieldOfView"]
acqFovCenter(f::MDFFileV1) = addLeadingSingleton(
              f["/acquisition/drivefield/fieldOfViewCenter"],2 )
acqFovCenter(f::MDFFileV2) = f["/acquisition/fieldOfViewCenter"]

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
# THE FOLLOWING IS NOT FULLY CORRECT!!!
rxBandwidth(f::MDFFileV1) =
    repeat( [f["/acquisition/receiver/bandwidth"]], outer=rxNumChannels(f))
rxBandwidth(f::MDFFileV2) = f["/acquisition/receiver/bandwidth"]
rxNumSamplingPoints(f::MDFFileV1) =
    repeat( [f["/acquisition/receiver/numSamplingPoints"]], outer=rxNumChannels(f))
rxNumSamplingPoints(f::MDFFileV2) = f["/acquisition/receiver/numSamplingPoints"]
function rxFrequencies(f::MDFFileV1)
  a = f["/acquisition/receiver/frequencies"]
  return reshape( repeat(a , outer=rxNumChannels(f)), length(a), rxNumChannels(f) )
end
rxFrequencies(f::MDFFileV2) = f["/acquisition/receiver/frequencies"]
rxTransferFunction(f::MDFFile) = f["/acquisition/receiver/transferFunction"]

# measurements
measUnit(f::MDFFileV1) = "a.u."
measUnit(f::MDFFileV2) = f["/measurement/unit"]
measRawDataConversion(f::MDFFileV1) = 1.0
measRawDataConversion(f::MDFFileV2) = f["/measurement/rawDataConversion"]
function measData(f::MDFFileV1)
  if !h5exists(f.filename, "/measurement")
    return nothing
  end
  tdExists = h5exists(f.filename, "/measurement/dataTD")

  if tdExists
    data = f["/measurement/dataTD"]
    if ndims(data) == 3
      return reshape(data,size(data,1),size(data,2),1,size(data,3))
    else
      return data
    end
  else
    data = f["/measurement/dataFD"]
    dataFD = reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4)))
    dataTD = irfft(dataFD, 2*(size(data,2)-1), 1)
    return reshape(dataTD,size(dataTD,1),size(dataTD,2),1,size(dataTD,3))
  end
end
measData(f::MDFFileV2) = f["/measurement/data"]
measDataTimeOrder(f::MDFFileV1) = collect(1:acqNumFrames(f))
measDataTimeOrder(f::MDFFileV2) = f["/measurement/dataTimeOrder"]
measBGData(f::MDFFile) = f["/measurement/backgroundData"]
measBGDataTimeOrder(f::MDFFile) = f["/measurement/backgroundDataTimeOrder"]

# calibrations
function calibSystemMatrixData(f::MDFFileV1)
  data = f["/calibration/dataFD"]
  if ndims(data) == 4
    return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),1))
  else
    return reinterpret(Complex{eltype(data)}, data, (size(data,2),size(data,3),size(data,4),size(data,5)))
  end
end
function calibSystemMatrixData(f::MDFFileV2)
  data = f["/calibration/systemMatrixData"]
  return reinterpret(Complex{eltype(data)}, data,
               (size(data,2),size(data,3),size(data,4),size(data,5)))
end
calibSNR(f::MDFFileV1) = f["/calibration/snrFD"]
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
recoData(f::MDFFileV1) = addTrailingSingleton(
         f[ "/reconstruction/data"], 3)
recoData(f::MDFFileV2) = f["/reconstruction/data"]
recoFov(f::MDFFile) = f["/reconstruction/fieldOfView"]
recoFovCenter(f::MDFFile) = f["/reconstruction/fieldOfViewCenter"]
recoSize(f::MDFFile) = f["/reconstruction/size"]
recoOrder(f::MDFFile) = f["/reconstruction/order"]
recoPositions(f::MDFFile) = f["/reconstruction/positions"]

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
