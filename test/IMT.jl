@testset "Testing IMT submodule" begin

# Download test files

fnMeas = "measurement3D.h5"
#fnSM = "SF3D.h5"

# Measurement File 

measIMT = MPIFile(fnMeas)
@test typeof(measIMT) == IMTFileMeas 

#calibIMT = MPIFile(fnSM)
#@test typeof(calibIMT) == IMTFileCalib

# only test this for v1
#@test uuid(mdfv1) == Base.Random.UUID("4b0ffb84-29f5-f388-49f2-92a206bba885")
#@test version(measIMT) == v"0.0.0"
#@test time(mdfv1) == DateTime("2016-02-08T14:28:34.673")
#@test experimentIsCalibration(calibIMT) == false
@test experimentIsCalibration(measIMT) == false
@test studyName(imt) == "n.a."

#for imt in (measIMT,calibIMT)
#  println("Test $imt")
#  @test studyName(imt) == "n.a."
#  @test studyNumber(imt) == 0
#  @test studyDescription(imt) == "n.a."

#  @test experimentName(imt) == "n.a."
#  @test experimentNumber(imt) == 0 
#  @test experimentDescription(imt) == "n.a."
#  @test experimentSubject(imt) == "n.a."
#  @test experimentIsSimulation(imt) == true 

#  @test scannerFacility(imt) == "n.a."
#  @test scannerOperator(imt) == "n.a."
#  @test scannerManufacturer(imt) == "n.a."
#  @test scannerName(imt) == "n.a."
#  @test scannerTopology(imt) == "n.a."

#  @test tracerName(imt) == ["n.a."]
#  @test tracerBatch(imt) == ["n.a."]
#  @test tracerVendor(imt) == ["n.a."]
#  @test tracerVolume(imt) == [0.0]
#  @test tracerConcentration(imt) == [0.0]
#  @test tracerInjectionTime(imt) == Dates.unix2datetime(0) #[DateTime("2015-09-15T11:17:23.011")]

  #@test acqStartTime(imt) == Dates.unix2datetime(0) #DateTime("2015-09-15T11:17:23.011")
  #@test acqGradient(imt)[:,:,1,1] == [-1.25 0 0; 0 -1.25 0;0 0 2.5]
  #@test acqFramePeriod(imt) == 6.528E-4
  #@test acqNumPeriodsPerFrame(imt) == 1
  #@test acqOffsetFieldShift(imt)[:,1,1] == [0.0; 0.0; -0.0]

  #@test dfNumChannels(imt) == 3
  #@test dfWaveform(imt) == "sine"
  #@test dfStrength(imt)[:,:,1] == [0.014 0.014 0.0]
  #@test dfPhase(imt)[:,:,1] == [1.5707963267948966 1.5707963267948966 1.5707963267948966]
  #@test dfBaseFrequency(imt) == 2500000.0
  #@test dfDivider(imt)[:,1] == [102; 96; 99]
  #@test dfCycle(imt) == 6.528E-4

  #@test rxNumChannels(imt) == 3
  #@test rxBandwidth(imt) == 1250000.0
  #@test rxNumSamplingPoints(imt) == 1632
  #@test acqNumAverages(imt) == 1

  #@test acqNumFrames(imt) == 500
  #@test acqNumPeriodsPerFrame(imt) == 1
  #@test acqNumPeriods(imt) == 500

  #@test size( measData(imt) ) == (1632,3,1,500)
  #@test size( measDataTDPeriods(imt) ) == (1632,3,500)
  #@test size( measDataTDPeriods(imt, 101:200) ) == (1632,3,100)

  #@test size(getMeasurements(imt, numAverages=1,
  #            spectralLeakageCorrection=false)) == (1632,3,1,500)

  #@test size(getMeasurements(imt, numAverages=10,
  #            spectralLeakageCorrection=false)) == (1632,3,1,50)

  #@test size(getMeasurements(imt, numAverages=10, frames=1:500,
  #            spectralLeakageCorrection=true)) == (1632,3,1,50)

  #@test size(getMeasurementsFD(imt, numAverages=10, frames=1:500)) == (817,3,1,50)

  #@test size(getMeasurementsFD(imt, numAverages=10, frames=1:500, loadasreal=true)) == (1634,3,1,50)

  #@test size(getMeasurementsFD(imt,frequencies=1:10, numAverages=10)) == (10,1,50)

#end



# Calibration File V1

#smv1 = MPIFile(fnSMV1)
#@test typeof(smv1) == MDFFileV1

#smv2 = MPIFile(fnSMV2)
#@test typeof(smv2) == MDFFileV2

#for sm in (smv1,smv2)
#  println("Test $sm")

#  @test size( measData(sm) ) == (1936,817,3,1)
#  @test measIsFourierTransformed(sm) == true
#  @test measIsTFCorrected(sm) == false
#  @test measIsTransposed(sm) == true
#  @test measIsBGCorrected(sm) == true

#  @test size( calibSNR(sm) ) == (817,3,1)
#  @test calibFov(sm) == [0.044; 0.044; 0.001]
#  @test calibFovCenter(sm) == [0.0; -0.0; 0.0]
#  @test calibSize(sm) == [44; 44; 1]
#  @test calibOrder(sm) == "xyz"
#  @test calibPositions(smv1) == nothing
#  @test calibOffsetField(smv1) == nothing
#  @test calibDeltaSampleSize(sm) == [0.001; 0.001; 0.001]
#  @test calibMethod(sm) == "robot"

#  @test size(filterFrequencies(sm, SNRThresh = 5)) == (147,)
#  #@test size(filterFrequencies(sm, numUsedFreqs = 100)) == (100,) # not working

#  @test size(getSystemMatrix(sm,1:10)) == (1936,10)
#  @test size(getSystemMatrix(sm,1:10,loadasreal=true)) == (1936,20)
#end

end
