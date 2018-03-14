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
  frames::Vector{Int64}
  offsets::Matrix{Float64}

  function MultiPatchMPIFile(filename::String, frames, offsets)
    return new(MPIFile(filename), frames, offsets)
  end

end


#function MultiPatchMPIFile(filename::String, roblog) #TODO
#
#end

function Base.show(io::IO, f::MultiPatchMPIFile)
  print(io, "Multi Patch MPI File: ", f.file)
end

acqNumPeriodsPerFrame(f::MultiPatchMPIFile) = length(f.frames)
acqNumFrames(f::MultiPatchMPIFile) = 1

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
  @eval $op(f::MultiPatchMPIFile) = $op(f.files[1])
end

for op in [ :dfStrength, :dfPhase ]
  @eval begin function $op(f::MultiPatchMPIFile)
       tmp = $op(f.file)
       newVal = similar(tmp, size(tmp,1), size(tmp,2),
                        acqNumPeriodsPerFrame(f))
       for y=1:acqNumPeriodsPerFrame(f)
         for a=1:size(tmp,1)
           for b=1:size(tmp,2)
             newVal[a,b,y] = tmp[a,b]
           end
         end
       end
      return newVal
    end
  end
end

function acqOffsetField(f::MultiPatchMPIFile)
   tmp = acqOffsetField(f.file)
   newVal = similar(tmp, 3, 1, acqNumPeriodsPerFrame(f))

   for b=1:acqNumPeriodsPerFrame(f)
     for a=1:3
         newVal[a,1,b] = tmp[a,1,1]
     end
   end

  return newVal
end

function acqGradient(f::MultiPatchMPIFile)
   tmp = acqGradient(f.file)
   newVal = similar(tmp, 3, 3, 1, acqNumPeriodsPerFrame(f))

   for b=1:acqNumPeriodsPerFrame(f)
     for a=1:3
       for d=1:3
           newVal[a,d,1,b] = tmp[a,d,1,1]
       end
     end
   end

  return newVal
end

for op in [:measIsFourierTransformed, :measIsTFCorrected,
           :measIsBGCorrected,
           :measIsTransposed, :measIsFramePermutation, :measIsFrequencySelection,
           :measIsSpectralLeakageCorrected,
           :measFramePermutation, :measIsBGFrame]
  @eval $op(f::MultiPatchMPIFile) = $op(f.files[1])
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
