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

acqNumPeriodsPerFrame(f::MultiMPIFile) = length(f.files)*acqNumFrames(f.files[1])
acqNumFrames(f::MultiMPIFile) = 1

for op in [:filepath, :version, :uuid, :time, :studyName, :studyNumber, :studyUuid, :studyDescription,
            :experimentName, :experimentNumber, :experimentUuid, :experimentDescription,
            :experimentSubject, :experimentHasMeasurement,
            :experimentIsSimulation, :experimentIsCalibration, :experimentHasProcessing,
            :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
            :tracerSolute, :tracerInjectionTime,
            :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName,
            :scannerTopology, :acqNumBGFrames, :acqFramePeriod,
            :acqStartTime,
            :dfNumChannels, :dfBaseFrequency, :dfDivider,
            :dfPeriod, :dfWaveform, :rxNumChannels, :acqNumAverages, :rxBandwidth,
            :rxNumSamplingPoints, :rxTransferFunction, :rxInductionFactor, :rxUnit, :rxDataConversionFactor]
  @eval $op(f::MultiMPIFile) = $op(f.files[1])
end

for op in [ :dfStrength, :dfPhase ]
  @eval begin function $op(f::MultiMPIFile)
       tmp = $op(f.files[1])
       newVal = similar(tmp, size(tmp,1), size(tmp,2), 
                        acqNumFrames(f.files[1]),length(f.files))
       for c=1:length(f.files)
         tmp = $op(f.files[c])
         for y=1:acqNumFrames(f.files[1])
           for a=1:size(tmp,1)
             for b=1:size(tmp,2)
               newVal[a,b,y,c] = tmp[a,b]
             end
           end
         end
       end
      return reshape(newVal,size(newVal,1),size(newVal,2),:)
    end
  end
end

for op in [ :acqGradient, :acqOffsetField, :acqOffsetFieldShift ]
  @eval begin function $op(f::MultiMPIFile)
       tmp = $op(f.files[1])
       newVal = similar(tmp, size(tmp,1), acqNumFrames(f.files[1]),length(f.files))
       for c=1:length(f.files)
         tmp = $op(f.files[c])
         for b=1:acqNumFrames(f.files[1])
           for a=1:size(tmp,1)
               newVal[a,b,c] = tmp[a]
           end
         end
       end
      return reshape(newVal,size(newVal,1),:)
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

##Achtung hack in der Schleife acqNumFrames(fi) statt acqNumFrames(f)
#notwendig, da hier Sprung zwischen MultiMPIFile und MPIFile
function measData(f::MultiMPIFile, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))
  data = zeros(Float64, rxNumSamplingPoints(f), length(receivers),
                        length(frames),length(periods),1)
  #for (i,p) in enumerate(periods)
  #  data[:,:,:,i,:] = measData(f.files[p], frames, 1, receivers)
  #end
  for (i,fi) in enumerate(f.files)
    fr_fi=acqNumFrames(fi)
    data[:,:,:,fr_fi*(i-1)+1:fr_fi*i,:] = measData(fi, 1:fr_fi, 1, receivers)
  end
  return reshape(data,size(data,1),size(data,2),:,1)
end



# TODO: define functions for multi calibration data
#  if experimentIsCalibration(f)
#    for op in [:calibSNR, :calibFov, :calibFovCenter,
#               :calibSize, :calibOrder, :calibPositions, :calibOffsetField,
#               :calibDeltaSampleSize, :calibMethod]
#      setparam!(params, string(op), eval(op)(f))
#    end
#  end
