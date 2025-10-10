@testset "getMeasurements Time and Frequency Domain" begin
  fnMeasBruker = joinpath(datadir,"BrukerStore","20150915_102110_Wuerfelphantom_1_1","18")
  fnSMBruker = joinpath(datadir,"BrukerStore","20141121_130749_CalibrationScans_1_1","76")


  # Preprocess SMs
  # Want to read from bruker, mdfv2 and mdfv2InMemory files
  # Data in time and frequency domain
  # With and without frequency filtering
  brukerMeas = MPIFile(fnMeasBruker)
  brukerSM = MPIFile(fnSMBruker)

  fnMPIMeas = joinpath(tmpdir, "getMeas", "meas.mdf")
  saveasMDF(fnMPIMeas, brukerMeas)
  mpiMeas = MPIFile(fnMPIMeas)
  memMeas = MDFv2InMemory(mpiMeas)

  fnMPISM = joinpath(tmpdir, "getMeas", "sm.mdf")
  saveasMDF(fnMPISM, brukerSM)
  mpiSM = MPIFile(fnMPISM)
  memSM = MDFv2InMemory(mpiSM)

  snr = 2
  fnFiltered = joinpath(tmpdir, "getMeas", "filtered.mdf")
  saveasMDF(fnFiltered, mpiSM; SNRThresh = snr)
  mpiFiltered = MPIFile(fnFiltered)
  memFiltered = MDFv2InMemory(mpiFiltered)
  freqs = vec([CartesianIndex(f,c) for f in measFrequencySelection(mpiFiltered), c in 1:rxNumChannels(mpiFiltered)])

  @testset "getMeasurements" begin
    # TODO
  end

  @testset "getMeasurementsFD" begin
    @testset "Measurement" begin
      # Default processing
      uBruker = getMeasurementsFD(brukerMeas)
      @test size(uBruker) == (817, 3, 1, acqNumFGFrames(brukerMeas))
      @test eltype(uBruker) == ComplexF32
      uMPI = getMeasurementsFD(mpiMeas)
      uMem = getMeasurementsFD(memMeas)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # With time domain processing such as spectralLeakageCorrection
      uBruker = getMeasurementsFD(brukerMeas, spectralLeakageCorrection = true)
      uMPI = getMeasurementsFD(mpiMeas, spectralLeakageCorrection = true)
      uMem = getMeasurementsFD(memMeas, spectralLeakageCorrection = true)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # With frequency filtering
      uBruker = getMeasurementsFD(brukerMeas, frequencies = freqs)
      @test size(uBruker) == (length(freqs), 1, acqNumFGFrames(brukerMeas))
      @test eltype(uBruker) == ComplexF32
      uMPI = getMeasurementsFD(mpiMeas, frequencies = freqs)
      uMem = getMeasurementsFD(memMeas, frequencies = freqs)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # Without time domain processing
      uAvgBruker = getMeasurementsFD(brukerMeas, frequencies = freqs, numAverages = acqNumFGFrames(brukerMeas))
      @test size(uAvgBruker) == (length(freqs), 1, 1)
      @test isapprox(mean(uBruker, dims = 3), uAvgBruker)
      uAvgMPI = getMeasurementsFD(mpiMeas, frequencies = freqs, numAverages = acqNumFGFrames(brukerMeas))
      uAvgMem = getMeasurementsFD(memMeas, frequencies = freqs, numAverages = acqNumFGFrames(brukerMeas))
      @test isapprox(uAvgBruker, uAvgMPI)
      @test isapprox(uAvgMPI, uAvgMem)
    end

    @testset "System Matrix" begin
      # Default processing
      uBruker = getMeasurementsFD(brukerSM, frames = 1:50)
      @test size(uBruker) == (817, 3, 1, 50)
      @test eltype(uBruker) == ComplexF32
      uMPI = getMeasurementsFD(mpiSM, frames = 1:50)
      uMem = getMeasurementsFD(memSM, frames = 1:50)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # With time domain processing such as spectralLeakageCorrection
      uBruker = getMeasurements(brukerSM, frames = 1:50, spectralLeakageCorrection = true)
      uMPI = getMeasurements(mpiSM, frames = 1:50, spectralLeakageCorrection = true)
      uMem = getMeasurements(memSM, frames = 1:50, spectralLeakageCorrection = true)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # With frequency filtering
      uBruker = getMeasurementsFD(brukerSM, frames = 1:50, frequencies = freqs)
      @test size(uBruker) == (length(freqs), 1, 50)
      @test eltype(uBruker) == ComplexF32
      uMPI = getMeasurementsFD(mpiSM, frames = 1:50, frequencies = freqs)
      uMem = getMeasurementsFD(memSM, frames = 1:50, frequencies = freqs)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)
      
      # Without time domain processing
      uAvgBruker = getMeasurementsFD(brukerSM, frames = 1:50, frequencies = freqs, numAverages = 50)
      @test size(uAvgBruker) == (length(freqs), 1, 1)
      @test isapprox(mean(uBruker, dims = 3), uAvgBruker)
      uAvgMPI = getMeasurementsFD(mpiSM, frames = 1:50, frequencies = freqs, numAverages = 50)
      uAvgMem = getMeasurementsFD(memSM, frames = 1:50, frequencies = freqs, numAverages = 50)
      @test isapprox(uAvgBruker, uAvgMPI)
      @test isapprox(uAvgMPI, uAvgMem)

      # With background frames
      bgFrames = measBGFrameIdx(brukerSM)
      @test_throws BoundsError getMeasurementsFD(brukerSM, frames = bgFrames, frequencies = freqs)
      uBruker = getMeasurementsFD(brukerSM, false, frames = bgFrames, frequencies = freqs)
      @test size(uBruker) == (length(freqs), 1, length(bgFrames))
      @test eltype(uBruker) == ComplexF32
      uMPI = getMeasurementsFD(mpiSM, false, frames = bgFrames, frequencies = freqs)
      uMem = getMeasurementsFD(memSM, false, frames = bgFrames, frequencies = freqs)
      @test isapprox(uBruker, uMPI)
      @test isapprox(uMPI, uMem)

      # Frequency filtered file
      uBruker = getMeasurementsFD(brukerSM, false, frames = bgFrames, frequencies = freqs)
      uMPIFiltered = getMeasurementsFD(mpiFiltered, false, frames = bgFrames)
      uMemFiltered = getMeasurementsFD(memFiltered, false, frames = bgFrames)
      @test size(uMPIFiltered) == (256, 3, 1, 23)
      @test isapprox(uBruker, reshape(uMPIFiltered, :, 1, 23))
      @test isapprox(reshape(uMPIFiltered, :, 1, 23), reshape(uMemFiltered, :, 1, 23))
      uMPIFiltered = getMeasurementsFD(mpiFiltered, false, frames = bgFrames, frequencies = freqs)
      uMemFiltered = getMeasurementsFD(memFiltered, false, frames = bgFrames, frequencies = freqs)
      @test isapprox(uBruker, uMPIFiltered)
      @test isapprox(uMPI, uMemFiltered)
      uBruker = getMeasurementsFD(brukerSM, false, frames = bgFrames, frequencies = freqs[1:2:end])
      uMPIFiltered = getMeasurementsFD(mpiFiltered, false, frames = bgFrames, frequencies = freqs[1:2:end])
      uMemFiltered = getMeasurementsFD(memFiltered, false, frames = bgFrames, frequencies = freqs[1:2:end])
      @test isapprox(uBruker, uMPIFiltered)
      @test isapprox(uMPIFiltered, uMemFiltered)
    end
  end
end