@testset "Frequency Filter" begin
  mdf = MDFv2InMemory()
  numChannels = 3
  bandwidth = 256
  numSamplingPoints = 2 * bandwidth
  recv = MDFv2Receiver(;bandwidth = bandwidth, numSamplingPoints = numSamplingPoints, numChannels = numChannels)
  mdf.acquisition = MDFv2Acquisition(;receiver = recv)
  mdf.measurement = MDFv2Measurement(;isFrequencySelection = false)

  unfiltered = filterFrequencies(mdf)
  @test length(unfiltered) == (bandwidth + 1) * numChannels
  @test eltype(unfiltered) == CartesianIndex{2}
  nFreq = bandwidth + 1

  @testset "Filter Receive Channels" begin
    @test length(filterFrequencies(mdf, recChannels = 1:numChannels)) == length(unfiltered)
    @test length(filterFrequencies(mdf, recChannels = 1:numChannels + 1)) == length(unfiltered)

    freqs = filterFrequencies(mdf, recChannels = 1:numChannels - 1)
    @test length(freqs) == length(unfiltered) - (nFreq)
    @test all(i -> in(i[2], 1:(numChannels - 1)), freqs)


    freqs = filterFrequencies(mdf, recChannels = 2:2)
    @test length(freqs) == nFreq
    @test all(i -> in(i[2], 2:2), freqs)

    chSelection = [1, 3]
    freqs = filterFrequencies(mdf, recChannels = chSelection)
    @test length(freqs) == length(unfiltered) - (nFreq)
    @test all(i -> in(i[2], chSelection), freqs)

    @test length(filterFrequencies(mdf, recChannels = 0:0)) == 0
    @test length(filterFrequencies(mdf, recChannels = numChannels+1:numChannels+2)) == 0
  end

  @testset "Min and Max Frequencies" begin
    # minFreq = 0 -> Offset Frequency with index 1
    @test length(filterFrequencies(mdf, minFreq = 0)) == length(unfiltered)
    @test length(filterFrequencies(mdf, minFreq = 1)) == length(unfiltered) - numChannels skip = true
    @test length(filterFrequencies(mdf, minFreq = bandwidth)) == numChannels skip = true
    @test length(filterFrequencies(mdf, minFreq = nFreq)) == 0

    @test length(filterFrequencies(mdf, maxFreq = bandwidth)) == length(unfiltered)
    @test length(filterFrequencies(mdf, maxFreq = bandwidth - 1)) == length(unfiltered) - numChannels skip = true
    @test length(filterFrequencies(mdf, maxFreq = 0)) == numChannels skip = true
  end

  @testset "Stop Bands" begin
    @test length(filterFrequencies(mdf, stopBands = [(0, 0)])) == length(unfiltered) - numChannels
    @test length(filterFrequencies(mdf, stopBands = [[0, 0]])) == length(unfiltered) - numChannels
    @test length(filterFrequencies(mdf, stopBands = [0:0])) == length(unfiltered) - numChannels

    @test length(filterFrequencies(mdf, stopBands = (0, 0))) == length(unfiltered) - numChannels
    @test length(filterFrequencies(mdf, stopBands = [0, 0])) == length(unfiltered) - numChannels
    @test length(filterFrequencies(mdf, stopBands = 0:0)) == length(unfiltered) - numChannels

    @test length(filterFrequencies(mdf, stopBands = 0:bandwidth)) == 0

    freqs = filterFrequencies(mdf, stopBands = 5:10)
    @test length(freqs) == length(unfiltered) - length(5:10) * numChannels
    @test !any(i -> in(i[1], (5:10) .+ 1), freqs)

    stopBands = [(5, 10), [15, 20], 18:22]
    freqs = filterFrequencies(mdf, stopBands = stopBands)
    @test length(freqs) == length(unfiltered) - (length(5:10) + length(15:22)) * numChannels
    @test !any(i -> in(i[1], (5:10) .+ 1) || in(i[1], (15:22) .+ 1), freqs)
  end


  snr = fill(0.0, nFreq, numChannels, 1)
  # Each freq has an SNR equal to its index
  snr[:, :, 1] .= 1:nFreq
  mdf.calibration = MDFv2Calibration(;snr = snr)

  @testset "SNR Threshold" begin
    @test length(filterFrequencies(mdf, SNRThresh = 0)) == length(unfiltered)
    @test length(filterFrequencies(mdf, SNRThresh = -1)) == length(unfiltered)
    @test length(filterFrequencies(mdf, SNRThresh = nFreq)) == numChannels
    @test length(filterFrequencies(mdf, SNRThresh = nFreq + 1)) == 0

    freqs = filterFrequencies(mdf, SNRThresh = floor(nFreq / 2))
    @test all(i -> i[1] >= floor(nFreq / 2), freqs)
  end

  @testset "Num Used Frequencies" begin
    @test length(filterFrequencies(mdf, numUsedFreqs = 1)) == numChannels skip = true
    @test length(filterFrequencies(mdf, numUsedFreqs = nFreq)) == length(unfiltered) skip = true
    @test length(filterFrequencies(mdf, numUsedFreqs = nFreq + 1)) == length(unfiltered) skip = true
    @test length(filterFrequencies(mdf, numUsedFreqs = nFreq - 1)) == length(unfiltered) - numChannels skip = true

    freqs = filterFrequencies(mdf, numUsedFreqs = 1)
    @test all(i -> i[1] == nFreq, freqs) skip = true
  end

end