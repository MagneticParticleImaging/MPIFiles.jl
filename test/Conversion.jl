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

  @testset "Frame Selection" begin
    frameSelectionMeas = joinpath(tmpdir, "conversion", "frames_meas.mdf")
    frames = 1:20:500
    saveasMDF(frameSelectionMeas, brukerMeas, frames = frames)
    MPIFile(frameSelectionMeas) do frameMeas
      @test isapprox(getMeasurements(brukerMeas, false, frames = frames), getMeasurements(frameMeas, false))
      @test size(measData(frameMeas), 4) == length(frames)
      @test length(measIsBGFrame(frameMeas)) == length(frames)
      @test measIsBGFrame(frameMeas) == measIsBGFrame(brukerMeas)[frames]
    end

    frameSelectionSM = joinpath(tmpdir, "conversion", "frames_sm.mdf")
    saveasMDF(frameSelectionSM, brukerSM, frames = frames)
    MPIFile(frameSelectionSM) do frameSM
      @test isapprox(getMeasurementsFD(brukerSM, false, frames = frames), getMeasurementsFD(frameSM, false))
      @test size(measData(frameSM), 1) == length(frames)
      @test length(measIsBGFrame(frameSM)) == length(frames)
      @test measIsBGFrame(frameSM) == measIsBGFrame(brukerSM)[frames]
      # Frame permutation should be in same "order" but refer to new frame subset
      @test sortperm(measFramePermutation(frameSM)) == sortperm(measFramePermutation(brukerSM)[frames])
      # TODO what should happen to calib fields?
    end
  end

  @testset "Frequency Filtering" begin
    @testset "Time Domain Origin" begin
    end
    @testset "Frequency Domain Origin" begin
    end
    @testset "Frequency Selected Origin" begin
    end
  end
end