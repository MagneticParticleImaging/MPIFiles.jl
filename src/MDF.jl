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

# general parameters
version(f::MDFFile) = VersionNumber( h5read(f.filename, "/version") )
uuid(f::MDFFile) = h5read(f.filename, "/uuid")
time(f::MDFFile) = DateTime( h5read(f.filename, "/date") )


# study parameters
studyName(f::MDFFile) = h5read(f.filename, "/study/name")
studyExperiment(f::MDFFileV1) = parse(Int64,h5read(f.filename, "/study/experiment"))
studyExperiment(f::MDFFileV2) = h5read(f.filename, "/study/experiment")
studyDescription(f::MDFFile) = h5read(f.filename, "/study/description")
studySubject(f::MDFFile) = h5read(f.filename, "/study/subject")
studyIsSimulation(f::MDFFileV2) = Bool( h5read(f.filename, "/study/isSimulation") )
studyIsSimulation(f::MDFFileV1) = Bool( h5read(f.filename, "/study/simulation") )
studyIsCalibration(f::MDFFileV2) = Bool( h5read(f.filename, "/study/isCalibration") )
function studyIsCalibration(f::MDFFileV1)
  return h5open(f.filename, "r") do file
    g = file["/"]
    exists(g, "calibration")
  end
end

# tracer parameters
tracerName(f::MDFFile) = h5read(f.filename, "/tracer/name")
tracerBatch(f::MDFFile) = h5read(f.filename, "/tracer/batch")
tracerVolume(f::MDFFile) = h5read(f.filename, "/tracer/volume")
tracerConcentration(f::MDFFile) = h5read(f.filename, "/tracer/concentration")
tracerSolute(f::MDFFileV2) = h5read(f.filename, "/tracer/solute")
tracerSolute(f::MDFFileV1) = "Fe"
tracerInjectionTime(f::MDFFileV1) = DateTime( h5read(f.filename, "/tracer/time") )
tracerInjectionTime(f::MDFFileV2) = DateTime( h5read(f.filename, "/tracer/injectionTime") )
tracerVendor(f::MDFFile) = h5read(f.filename, "/tracer/vendor")

# scanner parameters
scannerFacility(f::MDFFile) = h5read(f.filename, "/scanner/facility")
scannerOperator(f::MDFFile) = h5read(f.filename, "/scanner/operator")
scannerManufacturer(f::MDFFile) = h5read(f.filename, "/scanner/manufacturer")
scannerModel(f::MDFFile) = h5read(f.filename, "/scanner/model")
scannerTopology(f::MDFFile) = h5read(f.filename, "/scanner/topology")

# acquisition parameters
acqStartTime(f::MDFFileV1) = DateTime( h5read(f.filename, "/acquisition/time") )
acqStartTime(f::MDFFileV2) = DateTime( h5read(f.filename, "/acquisition/startTime") )
acqNumFrames(f::MDFFile) = h5read(f.filename, "/acquisition/numFrames")
acqNumBGFrames(f::MDFFileV2) = h5read(f.filename, "/acquisition/numBGFrames")
acqFramePeriod(f::MDFFile) = h5read(f.filename, "/acquisition/framePeriod")
acqNumPatches(f::MDFFile) = h5read(f.filename, "/acquisition/numPatches")
acqGradient(f::MDFFile) = addLeadingSingleton(
                              h5read(f.filename, "/acquisition/gradient"),2 )
acqOffsetField(f::MDFFileV2) = h5read(f.filename, "/acquisition/offsetField")
acqFov(f::MDFFileV1) = addLeadingSingleton(
               h5read(f.filename, "/acquisition/drivefield/fieldOfView"),2 )
acqFov(f::MDFFileV2) = h5read(f.filename, "/acquisition/fieldOfView")
acqFovCenter(f::MDFFileV1) = addLeadingSingleton(
              h5read(f.filename, "/acquisition/drivefield/fieldOfViewCenter"),2 )
acqFovCenter(f::MDFFileV2) = h5read(f.filename, "/acquisition/fieldOfViewCenter")

# drive-field parameters
dfNumChannels(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/numChannels")
dfStrength(f::MDFFileV1) = addTrailingSingleton( addLeadingSingleton(
         h5read(f.filename, "/acquisition/drivefield/strength"), 2), 3)
dfStrength(f::MDFFileV2) = h5read(f.filename, "/acquisition/drivefield/strength")
dfPhase(f::MDFFileV2) = h5read(f.filename, "/acquisition/drivefield/phase")
dfBaseFrequency(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/baseFrequency")
dfCustomWaveform(f::MDFFileV2) = h5read(f.filename, "/acquisition/drivefield/customWaveform")
dfDivider(f::MDFFileV1) = addTrailingSingleton(
                h5read(f.filename, "/acquisition/drivefield/divider"),2)
dfDivider(f::MDFFileV2) = h5read(f.filename, "/acquisition/drivefield/divider")
dfWaveform(f::MDFFileV2) = h5read(f.filename, "/acquisition/drivefield/waveform")
dfPeriod(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/period")

# receiver parameters
rxNumChannels(f::MDFFile) = h5read(f.filename, "/acquisition/receiver/numChannels")
rxNumAverages(f::MDFFileV1) = h5read(f.filename, "/acquisition/drivefield/averages")
rxNumAverages(f::MDFFileV2) = h5read(f.filename, "/acquisition/receiver/numAverages")
# THE FOLLOWING IS NOT FULLY CORRECT!!!
rxBandwidth(f::MDFFileV1) =
    repeat( [h5read(f.filename, "/acquisition/receiver/bandwidth")], outer=rxNumChannels(f))
rxBandwidth(f::MDFFileV2) = h5read(f.filename, "/acquisition/receiver/bandwidth")
rxNumSamplingPoints(f::MDFFileV1) =
    repeat( [h5read(f.filename, "/acquisition/receiver/numSamplingPoints")], outer=rxNumChannels(f))
rxNumSamplingPoints(f::MDFFileV2) = h5read(f.filename, "/acquisition/receiver/numSamplingPoints")
function rxFrequencies(f::MDFFileV1)
  a = h5read(f.filename, "/acquisition/receiver/frequencies")
  return reshape( repeat(a , outer=rxNumChannels(f)), length(a), rxNumChannels(f) )
end
rxFrequencies(f::MDFFileV2) = h5read(f.filename, "/acquisition/receiver/frequencies")
rxTransferFunction(f::MDFFile) = h5read(f.filename, "/acquisition/receiver/transferFunction")

# measurements
export measUnit, measRawDataConversion,
       measData, measDataTimeOrder, measBGData, measBGDataTimeOrder

# calibrations
export calibSystemMatrixData, calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

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
