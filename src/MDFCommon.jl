function measData(f::Union{MDFFileV2, MDFv2InMemory}, frames=1:acqNumFrames(f), periods=1:acqNumPeriodsPerFrame(f),
  receivers=1:rxNumChannels(f))

  _data = measDataRaw(f)
  if measIsFastFrameAxis(f)
    data = _data[frames, :, receivers, periods]
    data = reshape(data, length(frames), size(_data,2), length(receivers), length(periods))
  else
    data = _data[:, receivers, periods, frames]
    data = reshape(data, size(_data,1), length(receivers), length(periods), length(frames))
  end
  return data
end

function measDataTDPeriods(f::Union{MDFFileV2, MDFv2InMemory}, periods=1:acqNumPeriods(f),
  receivers=1:rxNumChannels(f))
  if measIsFastFrameAxis(f)
    error("measDataTDPeriods can currently not handle transposed data!")
  end

  _data = measDataRaw(f)
  data = reshape(_data,Val(3))[:, receivers, periods]

  return data
end

function systemMatrix(f::Union{MDFFileV2, MDFv2InMemory}, rows, bgCorrection=true)
  if !experimentHasMeasurement(f) || !measIsFourierTransformed(f)
    return nothing
  end

  rows_ = rowsToSubsampledRows(f, rows)

  data_ = measDataRaw(f)

  if !measIsFastFrameAxis(f)
    data_ = permutedims(data_, [4, 1, 2, 3])
  end

  data_ = data_[:, rows_, :]
  data = reshape(data_, Val(2))

  fgdata = data[measFGFrameIdx(f),:]

  if measIsSparsityTransformed(f)
    dataBackTrafo = similar(fgdata, prod(calibSize(f)), size(fgdata,2))

    B = createLinearOperator(measSparsityTransformation(f), ComplexF32; shape=tuple(calibSize(f)...))

    tmp = measSubsamplingIndices(f)
    subsamplingIndices_ = tmp[:, rows_, :]
    subsamplingIndices = reshape(subsamplingIndices_, Val(2))

    for l=1:size(fgdata,2)
      dataBackTrafo[:,l] .= 0.0
      dataBackTrafo[subsamplingIndices[:,l],l] .= fgdata[:,l]
      dataBackTrafo[:,l] .= adjoint(B) * vec(dataBackTrafo[:,l])
    end
    fgdata = dataBackTrafo
  end

  if bgCorrection # this assumes equidistant bg frames
    @debug "Applying bg correction on system matrix (MDF)"
    bgdata = data[measBGFrameIdx(f),:]
    blockLen = measBGFrameBlockLengths( invpermute!(deepcopy(measIsBGFrame(f)), measFramePermutation(f)) ) # Added deepcopy to be side-effect free in in-memory MDF
    st = 1
    for j=1:length(blockLen)
      bgdata[st:st+blockLen[j]-1,:] .=
      mean(bgdata[st:st+blockLen[j]-1,:], dims=1)
      st += blockLen[j]
    end

    bgdataInterp = interpolate(bgdata, (BSpline(Linear()), NoInterp()))
    # Cubic does not work for complex numbers
    origIndex = measFramePermutation(f)
    M = size(fgdata,1)
    K = size(bgdata,1)
    N = M + K
    for m=1:M
      alpha = (origIndex[m]-1)/(N-1)*(K-1)+1
      for k=1:size(fgdata,2)
        fgdata[m,k] -= bgdataInterp(alpha,k)
      end
    end
  end
  return fgdata
end

function systemMatrixWithBG(f::Union{MDFFileV2, MDFv2InMemory})
  if !experimentHasMeasurement(f) || !measIsFastFrameAxis(f) || !measIsFourierTransformed(f)
    return nothing
  end

  data_ = measDataRaw(f)
  data = data_[:, :, :, :]
  return data
end

# This is a special variant used for matrix compression
function systemMatrixWithBG(f::Union{MDFFileV2, MDFv2InMemory}, freq)
  if !experimentHasMeasurement(f) || !measIsFastFrameAxis(f) || !measIsFourierTransformed(f)
    return nothing
  end

  data_ = measDataRaw(f)
  data = data_[:, freq, :, :]
  return data
end