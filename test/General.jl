fnMeasBruker = "measurement_Bruker"
fnSMBruker = "systemMatrix_Bruker"
fnMeasV1 = "measurement_V1.mdf"
fnMeasV2 = "measurement_V2.mdf"
fnSMV1 = "systemMatrix_V1.mdf"
fnSMV2 = "systemMatrix_V2.mdf"

if !isdir(fnSMBruker)
  streamSM = get("http://media.tuhh.de/ibi/"*fnSMBruker*".zip")
  save(streamSM, fnSMBruker*".zip")
  run(`unzip $(fnSMBruker).zip`)
end
if !isdir(fnMeasBruker)
  streamMeas = get("http://media.tuhh.de/ibi/"*fnMeasBruker*".zip")
  save(streamMeas, fnMeasBruker*".zip")
  run(`unzip $(fnMeasBruker).zip`)
end


mdfBruker = MPIFile(fnMeasBruker)
@test typeof(mdfBruker) == BrukerFileMeas

saveasMDF(fnMeasV2, mdfBruker)#, frames=1:100) <- TODO test this

mdfv2 = MPIFile(fnMeasV2)
@test typeof(mdfv2) == MDFFileV2

for mdf in (mdfBruker,mdfv2)
  println("Test $mdf")
  @test studyName(mdf) == "Wuerfelphantom_Wuerfelphantom_1"
  @test studyNumber(mdf) == 1
  @test studyDescription(mdf) == "n.a."

  @test experimentName(mdf) == "fuenf (E18)"
  @test experimentNumber(mdf) == 18
  @test experimentDescription(mdf) == "fuenf (E18)"
  @test experimentSubject(mdf) == "Wuerfelphantom"
  @test experimentIsSimulation(mdf) == false
  @test experimentIsCalibration(mdf) == false

  @test scannerFacility(mdf) == "Universitätsklinikum Hamburg Eppendorf"
  @test scannerOperator(mdf) == "nmrsu"
  @test scannerManufacturer(mdf) == "Bruker/Philips"
  @test scannerModel(mdf) == "Preclinical MPI System"
  @test scannerTopology(mdf) == "FFP"

  @test tracerName(mdf) == ["Resovist"]
  @test tracerBatch(mdf) == ["0"]
  @test tracerVendor(mdf) == ["n.a."]
  @test tracerVolume(mdf) == [0.0]
  @test tracerConcentration(mdf) == [0.5]
  @test tracerInjectionTime(mdf) == [DateTime("2015-09-15T11:17:23.011")]

  @test acqStartTime(mdf) == DateTime("2015-09-15T11:17:23.011")
  @test acqGradient(mdf)[:,1] == [-1.25; -1.25; 2.5]
  @test acqFramePeriod(mdf) == 6.528E-4
  @test acqNumPatches(mdf) == 1
  @test acqOffsetFieldShift(mdf)[:,1] == [0.0; 0.0; -0.0]

  @test dfNumChannels(mdf) == 3
  @test dfWaveform(mdf) == "sine"
  @test dfStrength(mdf)[:,:,1] == [0.014 0.014 0.0]
  @test dfPhase(mdf)[:,:,1] == [1.5707963267948966 1.5707963267948966 1.5707963267948966]
  @test dfBaseFrequency(mdf) == 2500000.0
  @test dfDivider(mdf)[:,1] == [102; 96; 99]
  @test dfPeriod(mdf) == 6.528E-4

  @test rxNumChannels(mdf) == 3
  @test rxBandwidth(mdf) == 1250000.0
  @test rxNumSamplingPoints(mdf) == 1632
  @test rxNumAverages(mdf) == 1

  @test measNumFrames(mdf) == 500
  @test size( measData(mdf) ) == (1632,3,1,500)

  N = measNumFrames(mdf)

  @test size(getMeasurements(mdf, numAverages=1,
              spectralLeakageCorrection=false, fourierTransform=false)) == (1632,3,1,500)

  @test size(getMeasurements(mdf, numAverages=10,
              spectralLeakageCorrection=false, fourierTransform=false)) == (1632,3,1,50)

  @test size(getMeasurements(mdf, numAverages=10, frames=1:100,
              spectralLeakageCorrection=true, fourierTransform=false)) == (1632,3,1,10)

  @test size(getMeasurements(mdf, numAverages=10, frames=1:100,
              fourierTransform=true)) == (817,3,1,10)

  @test size(getMeasurements(mdf, numAverages=10, frames=1:100,
              fourierTransform=true, loadasreal=true)) == (1634,3,1,10)

  @test size(getMeasurements(mdf,frequencies=1:10, numAverages=10)) == (10,1,50)

end



# Calibration File

smBruker = MPIFile(fnSMBruker, isCalib=true)
@test typeof(smBruker) == BrukerFileCalib

saveasMDF(fnSMV2, smBruker)

smv2 = MPIFile(fnSMV2)
@test typeof(smv2) == MDFFileV2

for sm in (smBruker,smv2)
  println("Test $sm")

  @test size( systemMatrixWithBG(sm) ) == (1959,817,3,1)
  @test size( systemMatrix(sm,1:10) ) == (1936,10)

  # The following tests are commented out for the moment
  # since it is unclear what a BrukerFile should answer here
  #@test measIsFourierTransformed(sm) == true
  #@test measIsAveraged(sm) == false
  #@test measIsTFCorrected(sm) == false
  #@test measIsTransposed(sm) == true
  #@test measIsBGCorrected(sm) == true

  @test size( calibSNR(sm) ) == (817,3,1)
  @test calibFov(sm) == [0.044; 0.044; 0.001]
  @test calibFovCenter(sm) == [0.0; -0.0; 0.0]
  @test calibSize(sm) == [44; 44; 1]
  @test calibOrder(sm) == "xyz"
  @test calibDeltaSampleSize(sm) == nothing #[0.001; 0.001; 0.001]
  @test calibMethod(sm) == "robot"

  @test size(filterFrequencies(sm, SNRThresh = 5)) == (147,)
  #@test size(filterFrequencies(sm, numUsedFreqs = 100)) == (100,) # not working

  @test size(getSystemMatrix(sm,1:10)) == (1936,10)
  @test size(getSystemMatrix(sm,1:10,loadasreal=true)) == (1936,20)
end

# Next test checks if the cached system matrix is the same as the one loaded
# from the raw data
smBrukerPretendToBeMeas = MPIFile(fnSMBruker, isCalib=false)
S_loadedfromraw = getMeasurements(smBrukerPretendToBeMeas,
      frames=1:measNumFGFrames(smBrukerPretendToBeMeas),sortFrames=true,
      spectralLeakageCorrection=false,fourierTransform=true,transposed=true)

S_loadedfromproc = systemMatrix(smBruker)

@test norm(vec(S_loadedfromraw-S_loadedfromproc)) / norm(vec(S_loadedfromproc)) < 1e-6
