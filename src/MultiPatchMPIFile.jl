export MultiPatchMPIFile

#=
 TODO: - copy most of MultiMPIFile
       - The type should get
         1. the numberOfPeriodsPerPatch
         2. the periods at which these patches should be taken (initial periods
            of each patch)
         3. the new acqOffsetField (one has to take the gradient into account)

=#


type MultiPatchMPIFile <: MPIFile
  file::MPIFile
  periodRanges
  offsets

  function MultiPatchMPIFile(filename::String, periodRanges, offsets)
    return new(MPIFile(filename), periodRanges, offsets)
  end

end


#function MultiPatchMPIFile(filename::String, roblog) #TODO
#
#end

function Base.show(io::IO, f::MultiPatchMPIFile)
  print(io, "Multi Patch MPI File: ", f.file)
end

acqNumPeriodsPerFrame(f::MultiMPIFile) = length(f.periodRanges[1])
acqNumFrames(f::MultiMPIFile) = 1

for op in [:filepath, :version, :uuid, :time, :studyName, :studyNumber, :studyUuid, :studyDescription,
            :experimentName, :experimentNumber, :experimentUuid, :experimentDescription,
            :experimentSubject, :experimentHasMeasurement,
            :experimentIsSimulation, :experimentIsCalibration, :experimentHasProcessing,
            :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
            :tracerSolute, :tracerInjectionTime,
            :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName,
            :scannerTopology, :acqNumBGFrames,
            :acqStartTime,
            :dfNumChannels, :dfBaseFrequency, :dfDivider,
            :dfCycle, :dfWaveform, :rxNumChannels, :acqNumAverages, :rxBandwidth,
            :rxNumSamplingPoints, :rxTransferFunction, :rxInductionFactor, :rxUnit, :rxDataConversionFactor]
  @eval $op(f::MultiMPIFile) = $op(f.files[1])
end

for op in [ :dfStrength, :dfPhase ]
  @eval begin function $op(f::MultiPatchMPIFile)
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

function acqOffsetField(f::MultiPatchMPIFile)
   tmp = acqOffsetField(f.files[1])
   newVal = similar(tmp, 3, acqNumFrames(f.files[1]),length(f.files))
   for c=1:length(f.files)
     tmp = acqOffsetField(f.files[c])
     for b=1:acqNumFrames(f.files[1])
       for a=1:3
           newVal[a,b,c] = tmp[a,1,1,1]
       end
     end
   end
  return reshape(newVal,3,1,:)
end

function acqGradient(f::MultiPatchMPIFile)
   tmp = acqGradient(f.files[1])
   newVal = similar(tmp, 3, 3, acqNumFrames(f.files[1]),length(f.files))
   for c=1:length(f.files)
     tmp = acqGradient(f.files[c])
     for b=1:acqNumFrames(f.files[1])
       for a=1:3
         for d=1:3
             newVal[a,d,b,c] = tmp[a,d,1,1]
         end
       end
     end
   end
  return reshape(newVal,3,3,1,:)
end

for op in [:measIsFourierTransformed, :measIsTFCorrected,
           :measIsBGCorrected,
           :measIsTransposed, :measIsFramePermutation, :measIsFrequencySelection,
           :measIsSpectralLeakageCorrected,
           :measFramePermutation, :measIsBGFrame]
  @eval $op(f::MultiMPIFile) = $op(f.files[1])
end


experimentHasReconstruction(f::MultiPatchMPIFile) = false

##Achtung hack in der Schleife acqNumFrames(fi) statt acqNumFrames(f)
#notwendig, da hier Sprung zwischen MultiMPIFile und MPIFile

function measData(f::MultiPatchMPIFile, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))
  data = zeros(Float32, rxNumSamplingPoints(f), length(receivers),
                        length(frames),length(periods))
  #for (i,p) in enumerate(periods)
  #  data[:,:,:,i,:] = measData(f.files[p], frames, 1, receivers)
  #end
  for (i,fi) in enumerate(f.files)
    fr_fi=acqNumFrames(fi)
    data[:,:,:,fr_fi*(i-1)+1:fr_fi*i] = measData(fi, 1:fr_fi, 1, receivers)
  end
  return reshape(data,size(data,1),size(data,2),:,1)
end
