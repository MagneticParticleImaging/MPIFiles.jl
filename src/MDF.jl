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
date(f::MDFFile) = DateTime( h5read(f.filename, "/date") )


# study parameters
studyName(f::MDFFile) = h5read(f.filename, "/study/name")
studyExperiment(f::MDFFile) = h5read(f.filename, "/study/experiment")
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
tracerTime(f::MDFFile) = DateTime( h5read(f.filename, "/tracer/time") )

# scanner parameters
scannerFacility(f::MDFFile) = h5read(f.filename, "/scanner/facility")
scannerOperator(f::MDFFile) = h5read(f.filename, "/scanner/operator")
scannerManufacturer(f::MDFFile) = h5read(f.filename, "/scanner/manufacturer")
scannerModel(f::MDFFile) = h5read(f.filename, "/scanner/model")
scannerTopology(f::MDFFile) = h5read(f.filename, "/scanner/topology")

# acquisition parameters
export acqTime, acqNumFrames, acqNumBGFrames, acqFramePeriod, acqNumPatches,
       acqGradient, acqOffsetField, acqFov, acqFovCenter

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfPeriod

# receiver parameters
export rxNumChannels, rxNumAverages, rxBandwidth, rxNumSamplingPoints, rxFrequencies,
       rxTransferFunction

# measurements
export measUnit, measRawDataConversion,
       measData, measDataTimeOrder, measBGData, measBGDataTimeOrder

# calibrations
export calibSystemMatrixData, calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions





# study properties
subjectName(f::MDFFile) = h5read(f.filename, "/study/subject")
studyName(f::MDFFile) = h5read(f.filename, "/study/name")
measPath(f::MDFFile) = h5read(f.filename, "/study/_measPath")
expno(f::MDFFile) = h5read(f.filename, "/study/experiment")
description(f::MDFFile) = h5read(f.filename, "/study/description")
isReference(f::MDFFile) = Bool( h5read(f.filename, "/study/reference") )
isSimulation(f::MDFFile) = Bool( h5read(f.filename, "/study/simulation") )

# scanner properties

facility(f::MDFFile) = h5read(f.filename, "/scanner/facility")
operator(f::MDFFile) = h5read(f.filename, "/scanner/operator")
manufacturer(f::MDFFile) = h5read(f.filename, "/scanner/manufacturer")
model(f::MDFFile) = h5read(f.filename, "/scanner/model")
topology(f::MDFFile) = h5read(f.filename, "/scanner/topology")

# tracer properties
tracer(f::MDFFile) = h5read(f.filename, "/tracer/name")
tracerVendor(f::MDFFile)  = h5read(f.filename, "/tracer/vendor")
tracerBatch(f::MDFFile) = h5read(f.filename, "/tracer/batch")
tracerVolume(f::MDFFile) = h5read(f.filename, "/tracer/volume")
tracerConcentration(f::MDFFile) = h5read(f.filename, "/tracer/concentration")
timeTracerInjection(f::MDFFile) = DateTime( h5read(f.filename, "/tracer/time") )

# general acquisition properties
sfGradient(f::MDFFile) = h5read(f.filename, "/acquisition/gradient")
numScans(f::MDFFile) = h5read(f.filename, "/acquisition/numFrames")
acqDate(f::MDFFile) = DateTime( h5read(f.filename, "/acquisition/time") )
numPatches(f::MDFFile) = h5read(f.filename, "/acquisition/numPatches")
framePeriod(f::MDFFile) = h5read(f.filename, "/acquisition/framePeriod")

# drive-field properties
dfStrength(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/strength")
dfFov(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/fieldOfView")
dfcycle(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/period")
dfBaseFrequency(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/baseFrequency")
dfDivider(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/divider")
dfRepetitionTime(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/repetitionTime")
ffPos(f::MDFFile; alpha=[0,0,0]) = h5read(f.filename, "/acquisition/drivefield/fieldOfViewCenter") #FIXME
numAverages(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/averages")
numDfChannels(f::MDFFile) = h5read(f.filename, "/acquisition/drivefield/numChannels")

# receiver properties
bandwidth(f::MDFFile) = h5read(f.filename, "/acquisition/receiver/bandwidth")
frequencies(f::MDFFile) = h5read(f.filename, "/acquisition/receiver/frequencies")
numFreq(f::MDFFile) = div(numTimePoints(f),2)+1
numReceivers(f::MDFFile) = h5read(f.filename, "/acquisition/receiver/numChannels")
numTimePoints(f::MDFFile) = first(h5read(f.filename, "/acquisition/receiver/numSamplingPoints")) # komisch, Tobi fragen!

# system function properties
gridSize(f::MDFFile) = h5read(f.filename, "/calibration/size")
fov(f::MDFFile) = h5read(f.filename, "/calibration/fieldOfView")
deltaSampleConcentration(f::MDFFile) = 1.0

# measurements properties
samplePosition(f::MDFFile) = h5read(f.filename, "/measurement/samplePosition")
numSamplePosition(f::MDFFile) = h5read(f.filename, "/measurement/numSamplePosition")

filepath(f::MDFFile) = f.filename

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


export date
date(f::MDFFile) = h5read(f.filename, "/dates")
