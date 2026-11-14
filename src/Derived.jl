export acqNumFGFrames, acqNumBGFrames, acqOffsetFieldShift, acqFramePeriod,
       acqNumPeriods, acqNumPatches, acqNumPeriodsPerPatch, acqFov,
       acqGradientDiag,
       rxNumFrequencies, rxFrequencies, rxTimePoints,
       measFGFrameIdx, measBGFrameIdx, measBGFrameBlockLengths

"""
    rxNumFrequencies(f::MPIFile, numPeriodGrouping=1) -> Int

Number of frequencies of the real-valued Fourier transform of one drive-field period, or of
`numPeriodGrouping` concatenated periods.
"""
rxNumFrequencies(f::MPIFile, numPeriodGrouping=1) = length(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping))

"""
    rxFrequencies(f::MPIFile, numPeriodGrouping=1) -> Vector{Float64}

Frequencies in Hz corresponding to the frequency dimension of [`getMeasurementsFD`](@ref).
"""
rxFrequencies(f::MPIFile, numPeriodGrouping=1) = rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)) |> collect

"""
    rxTimePoints(f::MPIFile, numPeriodGrouping=1) -> Vector{Float64}

Sampling time points in s corresponding to the time dimension of [`getMeasurements`](@ref).
"""
rxTimePoints(f::MPIFile, numPeriodGrouping=1) = range(0, step=1/2rxBandwidth(f), length=rxNumSamplingPoints(f)*numPeriodGrouping) |> collect


"""
    acqGradientDiag(f::MPIFile) -> Array{Float64,3}

Diagonal entries of [`acqGradient`](@ref), size (3, Y, J).
"""
function acqGradientDiag(f::MPIFile)
  g = acqGradient(f)
  g_ = reshape(g,9,size(g,3),size(g,4))
  return g_[[1,5,9],:,:]
end

"""
    acqFov(f::MPIFile)

Field of view covered by the drive field in m, derived from [`dfStrength`](@ref) and [`acqGradient`](@ref).
"""
function acqFov(f::MPIFile)
  if size(dfStrength(f)[1,:,:],1) == 3
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[:,1,:] )
  else
    return  2*dfStrength(f)[1,:,:] ./ abs.( acqGradientDiag(f)[1,1,1] )
  end
end

"""
    acqFramePeriod(f::MPIFile) -> Float64

Duration of one frame in s, i.e. `dfCycle(f) * acqNumAverages(f) * acqNumPeriodsPerFrame(f)`.
"""
acqFramePeriod(b::MPIFile) = dfCycle(b) * acqNumAverages(b) * acqNumPeriodsPerFrame(b)

# numPeriods is the total number of DF periods in a measurement.
"""
    acqNumPeriods(f::MPIFile) -> Int

Total number of drive-field periods in the measurement, i.e. `acqNumFrames(f) * acqNumPeriodsPerFrame(f)`.
"""
acqNumPeriods(f::MPIFile) = acqNumFrames(f)*acqNumPeriodsPerFrame(f)

"""
    acqOffsetFieldShift(f::MPIFile) -> Array{Float64,3}

Spatial shift in m caused by the offset field, derived from [`acqOffsetField`](@ref) and [`acqGradient`](@ref).
"""
function acqOffsetFieldShift(f::MPIFile)
    return -acqOffsetField(f) ./ reshape( acqGradient(f),9,1,:)[[1,5,9],:,:]
end

"""
    acqNumFGFrames(f::MPIFile) -> Int

Number of foreground frames, denoted by O.
"""
acqNumFGFrames(f::MPIFile) = acqNumFrames(f) - acqNumBGFrames(f)
"""
    acqNumBGFrames(f::MPIFile) -> Int

Number of background frames, denoted by E.
"""
acqNumBGFrames(f::MPIFile) = sum(measIsBGFrame(f))

"""
    measBGFrameIdx(f::MPIFile) -> Vector{Int}

Indices of the background frames, see [`measIsBGFrame`](@ref).
"""
measBGFrameIdx(f::MPIFile) = findall(measIsBGFrame(f))
"""
    measFGFrameIdx(f::MPIFile) -> Vector{Int}

Indices of the foreground frames, see [`measIsBGFrame`](@ref).
"""
measFGFrameIdx(f::MPIFile) = findall(.!measIsBGFrame(f))

"""
    measBGFrameBlockLengths(mask) -> Vector{Int}

Lengths of the consecutive blocks of `true` values in the background frame `mask`.
"""
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


"""
    acqNumPatches(f::MPIFile) -> Int

Number of distinct patches, i.e. offset field positions. Not valid for varying gradients.
"""
function acqNumPatches(f::MPIFile)
  # not valid for varying gradients / multi gradient
  shifts = acqOffsetFieldShift(f)
  return size(unique(shifts,dims=3),3)
end

"""
    acqNumPeriodsPerPatch(f::MPIFile) -> Int

Number of drive-field periods per patch.
"""
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