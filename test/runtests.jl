using MPIFiles
using Base.Test
using Requests


# Download test files

fnMeasV1 = "measurement_V1.mdf"
fnSMV1 = "systemMatrix_V1.mdf"

if !isfile(fnSMV1)
  streamSM = get("http://media.tuhh.de/ibi/mdf/systemMatrix.h5")
  save(streamSM, fnSMV1)
end
if !isfile(fnMeasV1)
  streamMeas = get("http://media.tuhh.de/ibi/mdf/measurement_5.h5")
  save(streamMeas, fnMeasV1)
end

mdfv1 = MPIFile(fnMeasV1)
@test typeof(mdfv1) == MDFFileV1

@test uuid(mdfv1) == "4b0ffb8429f5f38849f292a206bba885"
@test version(mdfv1) == v"1.0.0"
@test time(mdfv1) == DateTime("2016-02-08T14:28:34.673")

@test studyName(mdfv1) == "Wuerfelphantom"
@test studyExperiment(mdfv1) == 18
@test studyDescription(mdfv1) == "n.a."
@test studySubject(mdfv1) == "Wuerfelphantom"
@test studyIsSimulation(mdfv1) == false
@test studyIsCalibration(mdfv1) == false

@test scannerFacility(mdfv1) == "University Medical Center Hamburg-Eppendorf, Germany"
@test scannerOperator(mdfv1) == "n.a."
@test scannerManufacturer(mdfv1) == "Bruker/Philips"
@test scannerModel(mdfv1) == "n.a."
@test scannerTopology(mdfv1) == "FFP"

@test tracerName(mdfv1) == "Resovist"
@test tracerBatch(mdfv1) == "0"
@test tracerVendor(mdfv1) == "n.a."
@test tracerVolume(mdfv1) == 0.0
@test tracerConcentration(mdfv1) == 0.5
@test tracerInjectionTime(mdfv1) == DateTime("2015-09-15T11:17:23.011")

@test acqStartTime(mdfv1) == DateTime("2015-09-15T11:17:23.011")
@test acqGradient(mdfv1) == [-1.25 -1.25 2.5]
@test acqFramePeriod(mdfv1) == 6.528E-4
@test acqNumFrames(mdfv1) == 500
@test acqNumPatches(mdfv1) == 1
@test acqFov(mdfv1) == [0.0224 0.0224 0.0]
@test acqFovCenter(mdfv1) == [0.0 0.0 -0.0]

@test dfNumChannels(mdfv1) == 3
@test dfStrength(mdfv1)[:,:,1] == [0.014 0.014 0.0]
@test dfBaseFrequency(mdfv1) == 2500000.0
@test dfDivider(mdfv1)[:,1] == [102; 96; 99]
@test dfPeriod(mdfv1) == 6.528E-4

@test rxNumChannels(mdfv1) == 3
@test rxBandwidth(mdfv1)[1] == 1250000.0
@test rxNumSamplingPoints(mdfv1)[1] == 1632
@test rxNumAverages(mdfv1) == 1
