export getMeasurements, getMeasurementsLowLevel

function measDataConv(f::MPIFile, args...)
  data = measData(f, args...)

  if eltype(data) <: Integer
    data = map(Float32, data)
  end
  a = rxDataConversionFactor(f)
  if a!=nothing
    if measIsTransposed(f)
      for d=1:size(data,3)
        slice = view(data,:,:,d,:)
        scale!(slice, a[1,d])
        slice .+= a[2,d]
      end
    else
      for d=1:size(data,2)
        slice = view(data,:,d,:,:)
        scale!(slice, a[1,d])
        slice .+= a[2,d]
      end
    end
  end
  return data
end

function spectralLeakageCorrection_(f::MPIFile, frames, periods)
  #println("Apply Spectral Cleaning")
  numTimePoints = rxNumSamplingPoints(f)
  numReceivers = rxNumChannels(f)
  numFrames = acqNumFrames(f)
  numPeriods = acqNumPeriods(f)

  data = zeros(Float32, numTimePoints, numReceivers, length(periods), length(frames))
  M = numTimePoints*3
  window3 = 0.5.*(1.-cos.(2*π/(M-1)*(0:(M-1))))
  window3 = window3 / (sum(window3)/M)

  M = numTimePoints*2
  window2 = 0.5.*(1.-cos.(2*π/(M-1)*(0:(M-1))))
  window2 = window2 / (sum(window2)/M)

  for (i,fr) in enumerate(frames)
    for (p,pe) in enumerate(periods)
      for r in 1:numReceivers
        if fr==1
          tmp = measDataConv(f, fr:fr+1, pe, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +  tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        elseif fr==numFrames
          tmp = measDataConv(f, fr-1:fr, pe, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +    tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        else
          tmp = measDataConv(f, fr-1:fr+1, pe, r)
          data[:,r,p,i] = 1/3 * (tmp[:,1,1,1] .* window3[1:numTimePoints]
                            +  tmp[:,1,1,2] .* window3[1+numTimePoints:2*numTimePoints]
                            +  tmp[:,1,1,3] .* window3[1+2*numTimePoints:3*numTimePoints]
                            );
        end
      end
    end
  end
  return data
end

function measDataLowLevel(f::MPIFile, args...; spectralLeakageCorrection=false )
  if measIsFourierTransformed(f)
    return measDataConv(f, args...)
  else
    if !spectralLeakageCorrection || measIsSpectralLeakageCorrected(f) ||
        acqNumFrames(f) == 1
       tmp = measDataConv(f, args...)
    else
       tmp = spectralLeakageCorrection_(f, args...)
    end
  end
end


function returnasreal{T}(u::AbstractArray{Complex{T}})
  return reinterpret(T,u,tuple(size(u,1)*2,size(u)[2:end]...))
end
returnasreal{T<:Real}(u::AbstractArray{T}) = u

function getAveragedMeasurements(f::MPIFile; frames=1:acqNumFrames(f),
            numAverages=1,  verbose = false, periods=1:acqNumPeriods(f),
            spectralLeakageCorrection=false)

  verbose && println( rxNumSamplingPoints(f), " ",
                      rxNumChannels(f), " ", acqNumFrames(f), )

  if numAverages == 1
    data = measDataLowLevel(f, frames, periods, spectralLeakageCorrection=spectralLeakageCorrection)
  else
    nFrames = length(frames)
    nBlocks = ceil(Int, nFrames / numAverages)

    rem(nFrames, numAverages) != 0 && (warn("numAverages no integer divisor of nFrames.
              Last Block will be averaged over less than $numAverages Frames."))

    if measIsTransposed(f)
      if measIsFourierTransformed(f)
        data = zeros(Complex64, nBlocks, rxNumFrequencies(f), rxNumChannels(f), acqNumPeriods(f))
      else
        data = zeros(Float32, nBlocks, rxNumSamplingPoints(f), rxNumChannels(f), acqNumPeriods(f))
      end
    else
      if measIsFourierTransformed(f)
        data = zeros(Complex64, rxNumFrequencies(f), rxNumChannels(f), acqNumPeriods(f), nBlocks)
      else
        data = zeros(Float32, rxNumSamplingPoints(f), rxNumChannels(f), acqNumPeriods(f), nBlocks)
      end
    end
    p = Progress(nBlocks, 1, "Loading measurement from $(filepath(f)) ...")
    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nFrames) # ensure that modulo is taken into account

      tmp = measDataLowLevel(f, frames[index1:index2], periods, spectralLeakageCorrection=spectralLeakageCorrection)
      if measIsTransposed(f)
        data[i,:,:,:] = mean(tmp,1)
      else
        data[:,:,:,i] = mean(tmp,4)
      end
      next!(p)
    end
  end
  return data
end

function getMeasurements(f::MPIFile, neglectBGFrames=true; frames=1:acqNumFrames(f),
      loadasreal=false, fourierTransform=measIsFourierTransformed(f),
      transposed=measIsTransposed(f), bgCorrection=false, frequencies=nothing,
      tfCorrection=measIsTFCorrected(f), sortFrames=false, kargs...)

  if neglectBGFrames
    idx = measFGFrameIdx(f)

    if !sortFrames
      data = getAveragedMeasurements(f; frames=idx[frames], kargs...)
    else
      data = getAveragedMeasurements(f; frames=idx[fgFramePermutation(f)][frames], kargs...)
    end

    if bgCorrection
      idxBG = measBGFrameIdx(f)
      dataBG = getAveragedMeasurements(f; frames=idxBG, kargs...)

      data[:,:,:,:] .-= mean(dataBG,4)
    end
  else
    if sortFrames
      perm1=cat(1,measFGFrameIdx(f),measBGFrameIdx(f))
      perm2=cat(1,fgFramePermutation(f),(length(perm1)-acqNumBGFrames(f)+1):length(perm1))
      permJoint = perm1[perm2]
      data = getAveragedMeasurements(f; frames=permJoint, kargs...)
    else
      data = getAveragedMeasurements(f; frames=frames, kargs...)
    end

    if bgCorrection
      idxBG = measBGFrameIdx(f)
      dataBG = getAveragedMeasurements(f; frames=idxBG, kargs...)

      data[:,:,:,:] .-= mean(dataBG,4)
    end
  end

  if (fourierTransform && !measIsFourierTransformed(f) ) || (frequencies != nothing)
    data = rfft(data, measIsTransposed(f) ? 2 : 1)
  end

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    if fourierTransform || measIsFourierTransformed(f) || (frequencies != nothing)
      data ./= tf
    else
      dim = measIsTransposed(f) ? 2 : 1
      J = size(data,dim)
      dataF = rfft(data, dim)
      dataF ./= tf
      data = irfft(dataF,J,dim)
    end
  end

  if frequencies != nothing
    # here we merge frequencies and channels
    if !measIsTransposed(f)
      data = reshape(data, size(data,1)*size(data,2), size(data,3), size(data,4))
      data = data[frequencies, :, :]
    else
      data = reshape(data, size(data,1), size(data,2)*size(data,3), size(data,4))
      data = data[:, frequencies, :]
    end
  end

  if transposed && !measIsTransposed(f)
    if frequencies != nothing
      data = permutedims(data, [3,1,2])
    else
      data = permutedims(data, [4,1,2,3])
    end
  end

  if loadasreal
    data = returnasreal(data)
  end

  return data
end
