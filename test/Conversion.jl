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
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 2
    params[:maxMixingOrder] = 5
    
    @testset "Time Domain Origin" begin
      freqSelectionTDOrigin = joinpath(tmpdir, "conversion", "freq_td.mdf")
      saveasMDF(joinpath(tmpdir, "conversion", "freq_td_1.mdf"), brukerMeas; params..., SNRThresh = -1) # measurements have no SNR
      MPIFile(freqSelectionTDOrigin) do freqTD
        @test measIsFrequencySelection(freqTD)
        @test measIsFourierTransformed(freqTD)
        freqs = [CartesianIndex(f,c) for f in measFrequencySelection(freqTD), c in 1:rxNumChannels(freqTD)]
        @test prod(size(measData(freqTD))[1:2]) == length(freqs)
        @test isapprox(getMeasurementsFD(brukerMeas, false, frequencies = freqs), getMeasurementsFD(freqTD, false))
      end
    end
    
    @testset "Frequency Domain Origin" begin
    end
    
    @testset "Frequency Selected Origin" begin
    end

    @testset "Frequency Selection Parameters & Handling" begin
    end
  end
end