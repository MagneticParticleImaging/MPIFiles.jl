@testset "Noise estimate" begin
  # Synthetic calibration data in the processed calibration layout
  # (frequency domain, fast frame axis, frame permutation), so that the
  # efficient code path of noiseEstimate is used.
  N, K, C, J = 12, 9, 2, 1
  isBG = fill(false, N)
  isBG[[1, 2, 6, 12]] .= true
  data = ComplexF32[sin(n + 2k + 3c) + im * cos(5n - k + c) for n = 1:N, k = 1:K, c = 1:C, j = 1:J]

  mdf = defaultMDFv2InMemory()
  rxNumChannels(mdf, C)
  rxNumSamplingPoints(mdf, 2(K - 1))
  rxBandwidth(mdf, 1.0e6)
  acqNumFrames(mdf, N)
  acqNumPeriodsPerFrame(mdf, J)
  measDataRaw(mdf, data)
  measIsFourierTransformed(mdf, true)
  measIsFastFrameAxis(mdf, true)
  measIsFramePermutation(mdf, true)
  measIsFrequencySelection(mdf, false)
  measIsBGFrame(mdf, isBG)
  measIsBGCorrected(mdf, false)
  measIsTFCorrected(mdf, false)
  measIsSpectralLeakageCorrected(mdf, false)
  measIsSparsityTransformed(mdf, false)
  @test measIsCalibProcessed(mdf)

  expected = std(data[isBG, :, :, 1], dims = 1)[1, :, :]
  freqs = [CartesianIndex(3, 1), CartesianIndex(5, 2), CartesianIndex(9, 1)]

  @testset "Background frames" begin
    n = noiseEstimate(mdf, tfCorrection = false)
    @test size(n) == (K, C)
    @test n ≈ expected
    @test noiseEstimate(mdf, frequencies = freqs, tfCorrection = false) ≈ expected[freqs]
  end

  @testset "Transfer function correction" begin
    tf = ComplexF64[1 + k / K + im * c / C for k = 1:K, c = 1:C]
    rxTransferFunction(mdf, tf)
    @test rxHasTransferFunction(mdf)
    # tfCorrection defaults to rxHasTransferFunction
    @test noiseEstimate(mdf) ≈ expected ./ abs.(tf)
    @test noiseEstimate(mdf, frequencies = freqs) ≈ expected[freqs] ./ abs.(tf[freqs])
    @test noiseEstimate(mdf, tfCorrection = false) ≈ expected

    inductionFactor = [2.0, -0.5]
    rxInductionFactor(mdf, inductionFactor)
    factor = abs.(tf) .* abs.(transpose(inductionFactor))
    @test noiseEstimate(mdf) ≈ expected ./ factor
    @test noiseEstimate(mdf, frequencies = freqs) ≈ expected[freqs] ./ factor[freqs]
    @test noiseEstimate(mdf, tfCorrection = false) ≈ expected

    # data that is already corrected is not corrected again
    measIsTFCorrected(mdf, true)
    @test noiseEstimate(mdf) ≈ expected
    measIsTFCorrected(mdf, false)
  end

  @testset "Stored noise estimate" begin
    fn = joinpath(tmpdir, "noiseEstimate.mdf")
    rxInductionFactor(mdf, nothing)
    saveasMDF(fn, mdf)
    stored = Float64[k + 10c for k = 1:K, c = 1:C]
    h5open(fn, "r+") do file
      file["/custom/noiseEstimate"] = stored
    end
    f = MPIFile(fn)
    @test noiseEstimate(f, tfCorrection = false) == stored
    @test noiseEstimate(f, frequencies = freqs, tfCorrection = false) == stored[freqs]
    # the stored estimate is uncorrected, the transfer function is applied when reading
    tf = rxTransferFunction(f)
    @test noiseEstimate(f) ≈ stored ./ abs.(tf)
    @test noiseEstimate(f, frequencies = freqs) ≈ stored[freqs] ./ abs.(tf[freqs])
    close(f)
  end
end
