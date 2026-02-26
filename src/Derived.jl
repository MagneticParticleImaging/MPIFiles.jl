export acqNumFGFrames, acqNumBGFrames, acqOffsetFieldShift, acqFramePeriod,
       acqNumPeriods, acqNumPatches, acqNumPeriodsPerPatch, acqFov,
       acqGradientDiag,
       rxNumFrequencies, rxFrequencies, rxTimePoints,
       measFGFrameIdx, measBGFrameIdx, measBGFrameBlockLengths, noiseEstimate

rxNumFrequencies(f::MPIFile, numPeriodGrouping=1) = length(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping))

rxFrequencies(f::MPIFile, numPeriodGrouping=1) = rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)) |> collect

rxTimePoints(f::MPIFile, numPeriodGrouping=1) = range(0, step=1/2rxBandwidth(f), length=rxNumSamplingPoints(f)*numPeriodGrouping) |> collect


function acqGradientDiag(f::MPIFile)
  g = acqGradient(f)
  g_ = reshape(g,9,size(g,3),size(g,4))
  return g_[[1,5,9],:,:]
end

function acqFov(f::MPIFile)
  if size(dfStrength(f)[1,:,:],1) == 3
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[:,1,:] )
  else
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[1,1,1] )
  end
end

acqFramePeriod(b::MPIFile) = dfCycle(b) * acqNumAverages(b) * acqNumPeriodsPerFrame(b)

# numPeriods is the total number of DF periods in a measurement.
acqNumPeriods(f::MPIFile) = acqNumFrames(f)*acqNumPeriodsPerFrame(f)

function acqOffsetFieldShift(f::MPIFile)
    return -acqOffsetField(f) ./ reshape( acqGradient(f),9,1,:)[[1,5,9],:,:]
end

acqNumFGFrames(f::MPIFile) = acqNumFrames(f) - acqNumBGFrames(f)
acqNumBGFrames(f::MPIFile) = sum(measIsBGFrame(f))

measBGFrameIdx(f::MPIFile) = findall(measIsBGFrame(f))
measFGFrameIdx(f::MPIFile) = findall(.!measIsBGFrame(f))

function measBGFrameBlockLengths(mask)
  len = Vector{Int}(undef,0)

  groupIdxStart = -1
  for i=1:(length(mask)+1)
    if i <= length(mask) && mask[i] && groupIdxStart == -1
      groupIdxStart = i
    end
    if groupIdxStart != -1 && ((i == length(mask)+1) || !mask[i])
      push!(len, i-groupIdxStart)
      groupIdxStart = - 1
    end
  end
  return len
end


function acqNumPatches(f::MPIFile)
  # not valid for varying gradients / multi gradient
  shifts = acqOffsetFieldShift(f)
  return size(unique(shifts,dims=3),3)
end

function acqNumPeriodsPerPatch(f::MPIFile)
  return div(acqNumPeriodsPerFrame(f), acqNumPatches(f))
end

export unflattenOffsetFieldShift

unflattenOffsetFieldShift(f::MPIFile) = unflattenOffsetFieldShift(acqOffsetFieldShift(f))
function unflattenOffsetFieldShift(shifts::Array)
  # not valid for varying gradients / multi gradient
  uniqueShifts = unique(shifts, dims=2)
  numPeriodsPerFrame = size(shifts,2)
  numUniquePatch = size(uniqueShifts,2)

  allPeriods = 1:numPeriodsPerFrame

  flatIndices = Vector{Vector{Int64}}()

  for i=1:numUniquePatch
    temp = allPeriods[vec(sum(shifts .== uniqueShifts[:,i],dims=1)).==3]
    push!(flatIndices, temp)
  end
  return flatIndices
end

# We assume that systemMatrixWithBG has already reordered the BG data
# to the end
systemMatrix(f::MPIFile) = systemMatrixWithBG(f)[1:acqNumFGFrames(f),:,:,:]

function measDataTD(f, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))

  data1 = measData(f,frames,periods,receivers)

  if measIsFastFrameAxis(f)
    data2 = permutedims(data1, invperm([4,1,2,3]))
  else
    data2 = data1
  end

  if measIsFourierTransformed(f)
    if measIsFrequencySelection(f)
      dataPadded = zeros(eltype(data2), (rxNumFrequencies(f), size(data2, 2), size(data2, 3), size(data2, 4)))
      dataPadded[measFrequencySelection(f), :, :, :] .= data2
      dataTD = irfft(dataPadded, rxNumSamplingPoints(f), 1)
    else
      dataTD = irfft(data2, rxNumSamplingPoints(f), 1)
    end
  else
    dataTD = data2
  end
  return dataTD
end
function measDataFD(f, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
                  receivers=1:rxNumChannels(f))
  data1 = measData(f,frames,periods,receivers)

  if measIsFastFrameAxis(f)
    data2 = permutedims(data1, invperm([4,1,2,3]))
  else
    data2 = data1
  end

  # TODO: frequencySelection
  if !measIsFourierTransformed(f)
    dataFD = rfft(data2, 1)
  else
    dataFD = data2
  end

  return dataFD
end

function noiseEstimate(f::MPIFile; frequencies=nothing, numPeriodGrouping=1, numPeriodAverages=1, kwargs...)
  if hasfield(typeof(f), :file) && haskey(f.file, "/custom/noiseEstimate") && numPeriodGrouping==1 && numPeriodAverages==1
    if !isnothing(frequencies)
      return f["/custom/noiseEstimate"][frequencies]
    else
      return f["/custom/noiseEstimate"]
    end
  else
    # Efficient shortcut useful for noise whitening during reconstruction
    if measIsCalibProcessed(f) && numPeriodGrouping==1 && numPeriodAverages==1
      data_ = measDataRaw(f)
      if !measIsFastFrameAxis(f)
        data_ = permutedims(data_, [4, 1, 2, 3])
      end
      if !isnothing(frequencies)
        return std(data_[measBGFrameIdx(f), frequencies, :],dims=1)[1,:,]
      else 
        return std(data_[measBGFrameIdx(f), :, :, :],dims=1)[1,:,:,]
      end
    end
    return std(getMeasurementsFD(f, false, frequencies=frequencies, frames=measBGFrameIdx(f), numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping, bgCorrection = false; kwargs...),dims=(3,4))[:,:,1,1]
  end
end