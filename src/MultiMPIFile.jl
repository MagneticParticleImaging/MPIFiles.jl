export MultiMPIFile

type MultiMPIFile <: MPIFile
  files::Vector{MPIFile}

  function MultiMPIFile(filenames::Vector{String})
    return new([MPIFile(f) for f in filenames])
  end

end

function Base.show(io::IO, f::MultiMPIFile)
  print(io, "Multi MPI File: ", f.files)
end

acqNumPeriods(f::MultiMPIFile) = length(f.files) #TODO make frames periods

for op in [:filepath, :version, :uuid, :time, :studyName, :studyNumber, :studyUuid, :studyDescription,
            :experimentName, :experimentNumber, :experimentUuid, :experimentDescription,
            :experimentSubject, :experimentHasMeasurement,
            :experimentIsSimulation, :experimentIsCalibration, :experimentHasProcessing,
            :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
            :tracerSolute, :tracerInjectionTime,
            :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName,
            :scannerTopology, :acqNumFrames, :acqNumBGFrames, :acqFramePeriod,
            :acqStartTime,
            :dfNumChannels, :dfBaseFrequency, :dfDivider,
            :dfPeriod, :dfWaveform, :rxNumChannels, :acqNumAverages, :rxBandwidth,
            :rxNumSamplingPoints, :rxTransferFunction, :rxInductionFactor, :rxUnit, :rxDataConversionFactor]
  @eval $op(f::MultiMPIFile) = $op(f.files[1])
end

for op in [ :dfStrength, :dfPhase ]
  @eval begin function $op(f::MultiMPIFile)
       tmp = $op(f.files[1])
       newVal = similar(tmp, size(tmp,1), size(tmp,2), acqNumPeriods(f))
       for c=1:acqNumPeriods(f)
         tmp = $op(f.files[c])
         for a=1:size(tmp,1)
           for b=1:size(tmp,2)
             newVal[a,b,c] = tmp[a,b]
           end
         end
       end
      return newVal
    end
  end
end

for op in [ :acqGradient, :acqOffsetField, :acqOffsetFieldShift ]
  @eval begin function $op(f::MultiMPIFile)
       tmp = $op(f.files[1])
       newVal = similar(tmp, size(tmp,1), acqNumPeriods(f))
       for c=1:acqNumPeriods(f)
         tmp = $op(f.files[c])
         for a=1:size(tmp,1)
             newVal[a,c] = tmp[a]
         end
       end
      return newVal
    end
  end
end

for op in [:measIsFourierTransformed, :measIsTFCorrected,
           :measIsBGCorrected,
           :measIsTransposed, :measIsFramePermutation, :measIsFrequencySelection,
           :measIsSpectralLeakageCorrected,
           :measFramePermutation, :measIsBGFrame]
  @eval $op(f::MultiMPIFile) = $op(f.files[1])
end


experimentHasReconstruction(f::MultiMPIFile) = false


function measData(f::MultiMPIFile, frames=1:acqNumFrames(f), periods=1:acqNumPeriods(f),
                  receivers=1:rxNumChannels(f))

  data = zeros(Float64, rxNumSamplingPoints(f), length(receivers),
                        length(periods), length(frames))
  for (i,p) in enumerate(periods)
    data[:,:,i,:] = measData(f.files[p], frames, 1, receivers)
  end
  return data
end



# TODO: define functions for multi calibration data
#  if experimentIsCalibration(f)
#    for op in [:calibSNR, :calibFov, :calibFovCenter,
#               :calibSize, :calibOrder, :calibPositions, :calibOffsetField,
#               :calibDeltaSampleSize, :calibMethod]
#      setparam!(params, string(op), eval(op)(f))
#    end
#  end
