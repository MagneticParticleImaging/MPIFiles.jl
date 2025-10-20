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
      saveasMDF(joinpath(tmpdir, "conversion", "freq_td.mdf"), brukerMeas; params..., SNRThresh = -1) # measurements have no SNR
      MPIFile(freqSelectionTDOrigin) do freqTD
        @test measIsFrequencySelection(freqTD)
        @test measIsFourierTransformed(freqTD)
        freqs = [CartesianIndex(f,c) for f in measFrequencySelection(freqTD), c in 1:rxNumChannels(freqTD)]
        @test prod(size(measData(freqTD))[1:2]) == length(freqs)
        @test isapprox(getMeasurementsFD(brukerMeas, false, frequencies = freqs), getMeasurementsFD(freqTD, false))
      end
    end
    
    @testset "Frequency Domain Origin" begin
      freqSelectionFDOrigin = joinpath(tmpdir, "conversion", "freq_fd.mdf")
      saveasMDF(joinpath(tmpdir, "conversion", "freq_fd.mdf"), brukerSM; params...)
      MPIFile(freqSelectionFDOrigin) do freqFD
        @test measIsFrequencySelection(freqFD)
        @test measIsFourierTransformed(freqFD)
        freqs = [CartesianIndex(f,c) for f in measFrequencySelection(freqFD), c in 1:rxNumChannels(freqFD)]
        @test prod(size(measData(freqFD))[2:3]) == length(freqs)
        @test isapprox(getMeasurementsFD(brukerSM, false, frequencies = freqs), getMeasurementsFD(freqFD, false))
      end
    end
    
    @testset "Frequency Selected Origin" begin
      freqSelectionTDOrigin = joinpath(tmpdir, "conversion", "freq_td.mdf")
      freqSelectionDoubleFiltered = joinpath(tmpdir, "conversion", "freq_double.mdf")
      if !isfile(freqSelectionTDOrigin)
        saveasMDF(joinpath(tmpdir, "conversion", "freq_td.mdf"), brukerMeas; params..., SNRThresh = -1) # measurements have no SNR
      end
      saveasMDF(freqSelectionDoubleFiltered, brukerMeas; maxMixingOrder = params[:maxMixingOrder] - 1)
      MPIFile(freqSelectionDoubleFiltered) do freqDouble
        @test measIsFrequencySelection(freqDouble)
        @test measIsFourierTransformed(freqDouble)
        freqs = [CartesianIndex(f,c) for f in measFrequencySelection(freqDouble), c in 1:rxNumChannels(freqDouble)]
        @test prod(size(measData(freqDouble))[1:2]) == length(freqs)
        @test isapprox(getMeasurementsFD(brukerMeas, false, frequencies = freqs), getMeasurementsFD(freqDouble, false))
      end
    end

    @testset "Frequency Selection Parameters & Handling" begin
      freqSelectionTDOrigin = joinpath(tmpdir, "conversion", "freq_td.mdf")
      freqSelectionDirectFiltered = joinpath(tmpdir, "conversion", "freq_direct.mdf")
      if !isfile(freqSelectionTDOrigin)
        saveasMDF(joinpath(tmpdir, "conversion", "freq_td.mdf"), brukerMeas; params..., SNRThresh = -1) # measurements have no SNR
      end
      freqs = [CartesianIndex(f,c) for f in measFrequencySelection(MPIFile(freqSelectionTDOrigin)), c in 1:rxNumChannels(brukerMeas)]
      
      
      @test_throws ArgumentError saveasMDF(freqSelectionDirectFiltered, MPIFile(freqSelectionTDOrigin); params..., frequencies = freqs)
      saveasMDF(freqSelectionDirectFiltered, MPIFile(freqSelectionTDOrigin); frequencies = freqs)
      MPIFile(freqSelectionDirectFiltered) do freqDirect
        @test measIsFrequencySelection(freqDirect)
        @test measIsFourierTransformed(freqDirect)
        @test prod(size(measData(freqDirect))[1:2]) == length(freqs)
        @test isapprox(getMeasurementsFD(MPIFile(freqSelectionTDOrigin), false), getMeasurementsFD(freqDirect, false))
      end

      # Check that apply frequency components to all channels equally
      freqs = [CartesianIndex(f,2) for f in measFrequencySelection(MPIFile(freqSelectionTDOrigin))]
      freqSelectionDirectExpanded = joinpath(tmpdir, "conversion", "freq_direct_expanded.mdf")
      saveasMDF(freqSelectionDirectExpanded, MPIFile(freqSelectionTDOrigin); frequencies = freqs)
      MPIFile(freqSelectionDirectExpanded) do freqExpanded
        @test isapprox(measData(freqExpanded), measData(MPIFile(freqSelectionTDOrigin)))
      end

      # Check that we can also only supply integer frequency components
      freqs = measFrequencySelection(MPIFile(freqSelectionTDOrigin))
      freqSelectionDirectIntegers = joinpath(tmpdir, "conversion", "freq_direct_integers.mdf")
      saveasMDF(freqSelectionDirectIntegers, MPIFile(freqSelectionTDOrigin); frequencies = freqs)
      MPIFile(freqSelectionDirectIntegers) do freqExpanded
        @test isapprox(measData(freqExpanded), measData(MPIFile(freqSelectionTDOrigin)))
      end
    end
  end
end