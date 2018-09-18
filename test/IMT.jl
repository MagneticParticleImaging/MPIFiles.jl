@testset "Testing IMT submodule" begin

# Download test files -> TODO

fnMeas = "measurement3D.h5"
fnCalib = "SF3D.h5"

# Measurement File 

measIMT = MPIFile(fnMeas)
calibIMT = MPIFile(fnCalib)
@test typeof(measIMT) == IMTFileMeas 
@test typeof(calibIMT) == IMTFileCalib 

#@test rxNumSamplingPoints(measIMT) == 53856 
@test rxNumSamplingPoints(calibIMT) == 26928 


for imt in (measIMT, measIMT)
  println("Test $imt")
  @test studyName(imt) == "n.a."
  @test studyNumber(imt) == 0
  @test studyDescription(imt) == "n.a."

  @test experimentName(imt) == "n.a."
  @test experimentNumber(imt) == 0 
  @test experimentDescription(imt) == "n.a."
  @test experimentSubject(imt) == "n.a."
  @test experimentIsSimulation(imt) == true 

  @test scannerFacility(imt) == "n.a." 
  @test scannerOperator(imt) == "n.a."
  @test scannerManufacturer(imt) == "n.a."
  @test scannerName(imt) == "n.a."
  @test scannerTopology(imt) == "n.a."

  #@test tracerName(imt) == ["n.a."]
  #@test tracerBatch(imt) == ["n.a."]
  @test tracerVendor(imt) == ["n.a."]
  @test tracerVolume(imt) == [0.0]
  @test tracerConcentration(imt) == [0.0]
  @test tracerInjectionTime(imt) == Dates.unix2datetime(0) 

  @test acqStartTime(imt) == Dates.unix2datetime(0) #DateTime("2015-09-15T11:17:23.011")
  @test acqGradient(imt)[:,:,1,1] == [0 0 0; 0 0 0; 0 0 0]
  #@test acqFramePeriod(imt) == 0.0215424
  @test acqNumPeriodsPerFrame(imt) == 1
  #@test acqOffsetFieldShift(imt)[:,1,1] == [0.0; 0.0; -0.0]
  @test acqNumAverages(imt) == 1
  @test acqNumFrames(imt) == 1
  #@test acqOffsetField(imt) == [0; 0; 0]
  #@test acqNumPeriods(imt) == 500

  @test dfNumChannels(imt) == 3
  @test dfWaveform(imt) == "sine"
  @test dfStrength(imt)[:,:,1] == [0.0 0.0 0.0]
  @test dfPhase(imt)[:,:,1] == [0.0 0.0 0.0] #[1.5707963267948966 1.5707963267948966 1.5707963267948966]
  @test dfBaseFrequency(imt) == 2500000.0
  @test dfDivider(imt)[:,1] == [102; 96; 99]
  @test dfCycle(imt) == [0.0215424]

  @test rxNumChannels(imt) == 3
  @test rxBandwidth(imt) == 1.25e6 
  #@test rxNumSamplingPoints(imt) == 53856 
  #@test rxUnit == "a.u."
  @test rxDataConversionFactor(imt) == reshape(Float64[1.0 0.0 1.0 0.0 1.0 0.0], 2,3) 

  @test size( measData(imt) ) == (53856,3,1,1)
  #@test size( measDataTDPeriods(imt) ) == (1632,3,500)
  #@test size( measDataTDPeriods(imt, 101:200) ) == (1632,3,100)

  #@test size(getMeasurements(imt, numAverages=1,
  #            spectralLeakageCorrection=false)) == (53856,3,1,1)

  #@test size(getMeasurements(imt, numAverages=10,
  #            spectralLeakageCorrection=false)) == (1632,3,1,50)

  #@test size(getMeasurements(imt, numAverages=10, frames=1:500,
  #            spectralLeakageCorrection=true)) == (1632,3,1,50)

  #@test size(getMeasurementsFD(imt, numAverages=10, frames=1:500)) == (817,3,1,50)

  #@test size(getMeasurementsFD(imt, numAverages=10, frames=1:500, loadasreal=true)) == (1634,3,1,50)

  #@test size(getMeasurementsFD(imt,frequencies=1:10, numAverages=10)) == (10,1,50)

end



# Calibration File V1

#for sm in (calibIMT)
#  println("Test $sm")

#  @test size( measData(sm) ) == (1936,817,3,1)
#  @test measIsFourierTransformed(sm) == true
#  @test measIsTFCorrected(sm) == false
#  @test measIsTransposed(sm) == true
#  @test measIsBGCorrected(sm) == false 

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
