
fnMeasBruker = "measurement_Bruker"
fnSMBruker = "systemMatrix_Bruker"


# Measurement File Bruker

b = MPIFile(fnMeasBruker)
@test typeof(b) == BrukerFile

@test studyName(b) == "1DSensitivity_1DSensitivity_1"
@test studyNumber(b) == 1
@test studyDescription(b) == "n.a."

@test experimentName(b) == "leer (E61)"
@test experimentNumber(b) == 61
@test experimentDescription(b) == "leer (E61)"
@test experimentSubject(b) == "1DSensitivity"
@test experimentIsSimulation(b) == false
@test experimentIsCalibration(b) == false

@test scannerFacility(b) == "Universitätsklinikum Hamburg Eppendorf"
@test scannerOperator(b) == "nmrsu"
@test scannerManufacturer(b) == "Bruker/Philips"
@test scannerModel(b) == "Preclinical MPI System"
@test scannerTopology(b) == "FFP"

@test tracerName(b)[1] == ""
@test tracerBatch(b)[1] == ""
@test tracerVendor(b)[1] == "n.a."
@test tracerVolume(b)[1] == 0.0
@test tracerConcentration(b)[1] == 0.5
@test tracerInjectionTime(b)[1] == DateTime("2015-01-20T10:52:31.897")

@test acqStartTime(b) == DateTime("2015-01-20T10:52:31.897")
@test acqGradient(b) == [-0.5 -0.5 1.0]
@test acqFramePeriod(b) == 2.15424
@test acqNumFrames(b) == 40
@test acqNumPatches(b) == 1
@test acqOffsetFieldShift(b) == [0.0 0.0 -0.0]

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

@test size(getMeasurements(b, numAverages=1,
            spectralLeakageCorrection=false, fourierTransform=false)) == (53856,3,1,40)

@test size(getMeasurements(b, numAverages=10,
            spectralLeakageCorrection=false, fourierTransform=false)) == (53856,3,1,4)

@test size(getMeasurements(b, numAverages=10,
            spectralLeakageCorrection=true, fourierTransform=false)) == (53856,3,1,4)

@test size(getMeasurements(b, numAverages=10,
            fourierTransform=true)) == (26929,3,1,4)

@test size(getMeasurements(b, numAverages=10,
            fourierTransform=true, loadasreal=true)) == (53858,3,1,4)


sm = MPIFile(fnSMBruker)
@test typeof(sm) == BrukerFile

@test experimentHasProcessing(sm) == true
@test size( procData(sm) ) == (100,26929,3,1)
@test procIsFourierTransformed(sm) == true
@test procIsAveraged(sm) == false
@test procIsTFCorrected(sm) == false
@test procIsTransposed(sm) == true
@test procIsBGCorrected(sm) == true

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
