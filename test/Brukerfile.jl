
fnMeasBruker = "measurement_Bruker"
fnSMBruker = "systemMatrix_Bruker"


# Measurement File Bruker

b = MPIFile(fnMeasBruker)
@test typeof(b) == BrukerFile

@test studyName(b) == "1DSensitivity_1DSensitivity_1"
@test studyExperiment(b) == 61
@test studyDescription(b) == "leer (E61)"
@test studySubject(b) == "1DSensitivity"
@test studyIsSimulation(b) == false
@test studyIsCalibration(b) == false

@test scannerFacility(b) == "Universitätsklinikum Hamburg Eppendorf"
@test scannerOperator(b) == "nmrsu"
@test scannerManufacturer(b) == "Bruker/Philips"
@test scannerModel(b) == "Preclinical MPI System"
@test scannerTopology(b) == "FFP"

@test tracerName(b) == ""
@test tracerBatch(b) == ""
@test tracerVendor(b) == "n.a."
@test tracerVolume(b) == 0.0
@test tracerConcentration(b) == 0.5
@test tracerInjectionTime(b) == DateTime("2015-01-20T10:52:31.897")

@test acqStartTime(b) == DateTime("2015-01-20T10:52:31.897")
@test acqGradient(b) == [-0.5 -0.5 1.0]
@test acqFramePeriod(b) == 2.15424
@test acqNumFrames(b) == 40
@test acqNumPatches(b) == 1
@test acqFov(b) == [0.0 0.0 0.028]
@test acqFovCenter(b) == [0.0 0.0 -0.0]

@test dfNumChannels(b) == 3
@test dfWaveform(b) == "sine"
@test dfStrength(b)[:,:,1] == [0.0 0.0 0.014]
@test dfPhase(b)[:,:,1] == [1.5707963267948966 1.5707963267948966 1.5707963267948966]
@test dfBaseFrequency(b) == 2500000.0
@test dfDivider(b)[:,1] == [102; 96; 99]
@test dfPeriod(b) == 0.0215424

@test rxNumChannels(b) == 3
@test rxBandwidth(b) == 1250000.0
@test rxNumSamplingPoints(b) == 53856
@test rxNumAverages(b) == 100

@test size( measData(b) ) == (53856,3,1,40)


sm = MPIFile(fnSMBruker)
@test typeof(sm) == BrukerFile

@test size( calibSystemMatrixData(sm) ) == (100,26929,3,1)
@test size( calibSNR(sm) ) == (26929,3)
@test calibFov(sm) == [0.001,0.001,0.04]
@test calibFovCenter(sm) == [0.0; -0.0; 0.0]
@test calibSize(sm) == [1,1,100]
@test calibOrder(sm) == "xyz"
@test calibPositions(sm) == nothing
@test calibOffsetField(sm) == nothing
@test calibDeltaSampleSize(sm) == nothing #TODO
@test calibMethod(sm) == "robot"

fnMeasConv = "measurement_conv.mdf"
fnSMConv = "systemMatrix_conv.mdf"

saveasMDF(fnMeasConv, fnMeasBruker)
saveasMDF(fnSMConv, fnSMBruker)
