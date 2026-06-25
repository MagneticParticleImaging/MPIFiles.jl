using Dates
using HTTP
using HDF5
using LinearAlgebra
using MPIFiles
using Statistics
using Test
using UUIDs
using Unitful
using Scratch
using LazyArtifacts
using Aqua
using Dagger

const datadir = joinpath(artifact"data", "data")
@info "The test data is located at $datadir."

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."

function compareMetadata(a::MPIFile, b::MPIFile; calib = true, meas = true, acq = true, df = true, rx = true, tracer = true, scanner = true, reco = true)
  if scanner
    @test scannerFacility(a) == scannerFacility(b)
    @test scannerOperator(a) == scannerOperator(b)
    @test scannerManufacturer(a) == scannerManufacturer(b)
    @test scannerName(a) == scannerName(b)
    @test scannerTopology(a) == scannerTopology(b)
  end

  if tracer
    @test tracerName(a) == tracerName(b)
    @test tracerBatch(a) == tracerBatch(b)
    @test tracerVendor(a) == tracerVendor(b)
    @test tracerVolume(a) == tracerVolume(b)
    @test tracerConcentration(a) == tracerConcentration(b)
    @test tracerInjectionTime(a) == tracerInjectionTime(b)
  end

  if acq
    @test acqStartTime(a) == acqStartTime(b)
    @test acqGradient(a) == acqGradient(b)
    @test acqFramePeriod(a) == acqFramePeriod(b)
    @test acqNumPeriodsPerFrame(a) == acqNumPeriodsPerFrame(b)
    @test acqOffsetFieldShift(a) == acqOffsetFieldShift(b)
  end

  if df
    @test dfNumChannels(a) == dfNumChannels(b)
    @test dfWaveform(a) == dfWaveform(b)
    @test dfStrength(a) == dfStrength(b)
    @test dfPhase(a) == dfPhase(b)
    @test dfBaseFrequency(a) == dfBaseFrequency(b)
    @test dfDivider(a) == dfDivider(b)
    @test dfDivider(a) == dfDivider(b)
  end

  if rx
    @test rxNumChannels(a) == rxNumChannels(b)
    @test rxBandwidth(a) == rxBandwidth(b)
    @test rxNumSamplingPoints(a) == rxNumSamplingPoints(b)
    @test rxDataConversionFactor(a) == rxDataConversionFactor(b)
    @test rxTransferFunction(a) == rxTransferFunction(b)
    @test rxTransferFunctionFileName(a) == rxTransferFunctionFileName(b)
  end

  if meas
    @test measIsFourierTransformed(a) == measIsFourierTransformed(b)
    @test measIsTFCorrected(a) == measIsTFCorrected(b)
    @test measIsBGCorrected(a) == measIsBGCorrected(b)
    @test rxDataConversionFactor(a) == rxDataConversionFactor(b)
    @test measIsFastFrameAxis(a) == measIsFastFrameAxis(b)
    @test measIsFramePermutation(a) == measIsFramePermutation(b)
    @test measIsFrequencySelection(a) == measIsFrequencySelection(b)
    @test measIsSpectralLeakageCorrected(a) == measIsSpectralLeakageCorrected(b)
    @test measFramePermutation(a) == measFramePermutation(b)
    @test measIsBGFrame(a) == measIsBGFrame(b)
    @test measIsFrequencySelection(a) == measIsFrequencySelection(b)
  end

  if calib
    @test calibSNR(a) == calibSNR(b)
    @test calibFov(a) == calibFov(b)
    @test calibFovCenter(a) == calibFovCenter(b)
    @test calibSize(a) == calibSize(b)
    @test calibOrder(a) == calibOrder(b)
    @test calibDeltaSampleSize(a) == calibDeltaSampleSize(b)
    @test calibMethod(a) == calibMethod(b)
    @test calibPositions(a) == calibPositions(b)
  end

  if reco
    @test recoFieldOfView(a) == recoFieldOfView(b)
    @test recoFieldOfViewCenter(a) == recoFieldOfViewCenter(b)
    @test recoPositions(a) == recoPositions(b)
    @test recoOrder(a) == recoOrder(b)
    @test recoSize(a) == recoSize(b)
  end

end

mkpath(joinpath(tmpdir,"mdf"))
mkpath(joinpath(tmpdir,"mdfim"))
mkpath(joinpath(tmpdir, "getMeas"))
mkpath(joinpath(tmpdir,"positions"))
mkpath(joinpath(tmpdir,"transferFunction"))
mkpath(joinpath(tmpdir,"conversion"))

@testset "Aqua" begin
  Aqua.test_all(MPIFiles)
end

include("General.jl")
include("Conversion.jl")
include("Cartesian.jl")
include("Positions.jl")
include("Subsampling.jl")
include("Interpolation.jl")
include("MDFv1.jl")
include("MultiMPIFile.jl")
include("Reco.jl")
include("IMT.jl")
include("CustomSFMeas.jl")
include("MDFInMemory.jl")
include("TransferFunction.jl")
include("FrequencyFilter.jl")
include("MagneticFieldMeasurement.jl")
include("DatasetStore.jl")
include("Measurements.jl")

@info "The unit tests are done!"
