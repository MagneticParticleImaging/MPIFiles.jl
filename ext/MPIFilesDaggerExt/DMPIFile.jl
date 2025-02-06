struct DaggerMPIFile{T, C <: Dagger.Chunk{T}} <: DMPIFile{T}
  file::C
  worker::Int64
end
function MPIFiles.DMPIFile(args...; worker::Int64, kwargs...)
  chunk = Dagger.@mutable worker = worker MPIFile(args...; kwargs...)
  return DaggerMPIFile(chunk, worker)
end  
function MPIFiles.DMPIFile(mdf::MDFv2InMemory; worker::Int64)
  chunk = Dagger.@mutable worker = worker mdf
  return DaggerMPIFile(chunk, worker)
end
MPIFiles.worker(dmdf::DaggerMPIFile) = dmdf.worker


field_symbols = [
# general
:version, :uuid,
# study parameters
:studyName, :studyNumber, :studyUuid, :studyDescription, :studyTime,
# experiment parameters
:experimentName, :experimentNumber, :experimentUuid, :experimentDescription,
:experimentSubject, :experimentIsSimulation, :experimentIsCalibration,
:experimentHasMeasurement, :experimentHasReconstruction,
# tracer parameters
:tracerName, :tracerBatch, :tracerVolume, :tracerConcentration, :tracerSolute,
:tracerInjectionTime, :tracerVendor,
# scanner parameters
:scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName, :scannerTopology,
# acquisition parameters
:acqStartTime, :acqNumFrames, :acqNumAverages, :acqGradient, :acqOffsetField,
:acqNumPeriodsPerFrame,
# drive-field parameters
:dfNumChannels, :dfStrength, :dfPhase, :dfBaseFrequency, :dfCustomWaveform, :dfDivider,
:dfWaveform, :dfCycle,
# receiver parameters
:rxNumChannels, :rxBandwidth, :rxNumSamplingPoints, :rxTransferFunction, :rxUnit,
:rxDataConversionFactor, :rxInductionFactor, :rxHasTransferFunction, :rxTransferFunctionFileName,
# measurements
:measData, :measDataTDPeriods, :measIsFourierTransformed, :measIsTFCorrected,
:measIsBGCorrected, :measIsFastFrameAxis, :measIsFramePermutation, :measIsFrequencySelection,
:measIsBGFrame, :measIsSpectralLeakageCorrected, :measFramePermutation,
# calibrations
:calibSNR, :calibFov, :calibFovCenter, :calibSize, :calibOrder, :calibPositions,
:calibOffsetFields, :calibDeltaSampleSize, :calibMethod, :calibIsMeanderingGrid,
# reconstruction results
:recoData, :recoFov, :recoFovCenter, :recoSize, :recoOrder, :recoPositions,
# additional functions that should be implemented by an MPIFile
:filepath, :systemMatrixWithBG, :systemMatrix, :selectedChannels,
]

for field in field_symbols
  @eval begin
    function MPIFiles.$field(file::DaggerMPIFile, args...; kwargs...)
      fetch(Dagger.spawn(file.file) do f
        $field(f, args...; kwargs...)
      end)
    end
  end
end

for fun in [:getMeasurements, :getMeasurementsFD, :filterFrequencies, :sortFrequencies]
  @eval begin
    function MPIFiles.$fun(file::DaggerMPIFile, args...; kwargs...)
      fetch(Dagger.spawn(file.file) do f
        $fun(f, args...; kwargs...)
      end)
    end
  end
end