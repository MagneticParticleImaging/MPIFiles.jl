@testset "saveasMDF" begin
  fnMeasBruker = joinpath(datadir,"BrukerStore","20150915_102110_Wuerfelphantom_1_1","18")
  fnSMBruker = joinpath(datadir,"BrukerStore","20141121_130749_CalibrationScans_1_1","76")

  brukerMeas = MPIFile(fnMeasBruker)
  brukerSM = MPIFile(fnSMBruker)

  fnMPIMeas = joinpath(tmpdir, "conversion", "meas.mdf")
  saveasMDF(fnMPIMeas, brukerMeas)
  mpiMeas = MPIFile(fnMPIMeas)
  compareMetadata(brukerMeas, mpiMeas; reco = false, calib = false)
  @test isapprox(measData(brukerMeas), measData(mpiMeas))

  fnMPISM = joinpath(tmpdir, "conversion", "sm.mdf")
  saveasMDF(fnMPISM, brukerSM)
  mpiSM = MPIFile(fnMPISM)
  compareMetadata(brukerSM, mpiSM; reco = false)
  @test isapprox(measData(brukerSM), measData(mpiSM))


  @testset "Frequency Filtering" begin
    @testset "Time Domain Origin" begin
    end
    @testset "Frequency Domain Origin" begin
    end
    @testset "Frequency Selected Origin" begin
    end
  end
end