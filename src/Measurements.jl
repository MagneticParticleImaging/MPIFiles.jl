export getMeasurements, getMeasurementsLowLevel, getMeasurementsFT

function spectralLeakageCorrection_(f::MPIFile, frames)
  #println("Apply Spectral Cleaning")
  numTimePoints = rxNumSamplingPoints(f)
  numReceivers = rxNumChannels(f)
  numFrames = acqNumFrames(f)
  numPatches = acqNumPatches(f)

  data = zeros(numTimePoints, numReceivers, numPatches, length(frames))
  M = numTimePoints*3
  window3 = 0.5.*(1.-cos(2*π/(M-1)*(0:(M-1))))
  window3 = window3 / (sum(window3)/M)

  M = numTimePoints*2
  window2 = 0.5.*(1.-cos(2*π/(M-1)*(0:(M-1))))
  window2 = window2 / (sum(window2)/M)

  for (i,fr) in enumerate(frames)
    for p in 1:numPatches
      for r in 1:numReceivers
        if fr==1
          tmp = measDataConv(f, fr:fr+1, p, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +  tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        elseif fr==numFrames
          tmp = measDataConv(f, fr-1:fr, p, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +    tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        else
          tmp = measDataConv(f, fr-1:fr+1, p, r)
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

function returnasreal{T}(u::AbstractArray{Complex{T}})
  return reinterpret(T,u,tuple(size(u,1)*2,size(u)[2:end]...))
end
returnasreal{T<:Real}(u::AbstractArray{T}) = u

function getMeasurementsLowLevel(f::MPIFile; frames=1:acqNumAllFrames(f),
            numAverages=1,  verbose = false,
            spectralLeakageCorrection=true)

  verbose && println( rxNumSamplingPoints(f), " ",
                      rxNumChannels(f), " ", acqNumFrames(f), )

  if numAverages == 1
    if !spectralLeakageCorrection # || hasSpectralCleaning(f)
      data = measDataConv(f, frames)
    else
      data = spectralLeakageCorrection_(f, frames)
    end
  else
    nFrames = length(frames)
    nBlocks = ceil(Int, nFrames / numAverages)

    rem(nFrames, numAverages) != 0 && (warn("numAverages no integer divisor of nFrames.
              Last Block will be averaged over less than $numAverages Frames."))

    data = zeros(rxNumSamplingPoints(f),rxNumChannels(f), acqNumPatches(f), nBlocks)

    p = Progress(nBlocks, 1, "Loading measurement from $(filepath(f)) ...")
    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nFrames) # ensure that modulo is taken into account

      if !spectralLeakageCorrection #|| hasSpectralCleaning(f)
        tmp = measDataConv(f, frames[index1:index2])
      else
        tmp = spectralLeakageCorrection_(f, frames[index1:index2])
      end
      data[:,:,:,i] = mean(tmp,4)
      next!(p)
    end
  end
  return data
end

function getMeasurementsLowLevelFromProc(f::MPIFile; frames=1:acqNumFrames(f),
            numAverages=1, spectralLeakageCorrection=true)

  data = procData(f, frames)
  if procIsFourierTransformed(f)
    data = irfft(data, 2*(size(data,1)-1), 1)
  end

  if numAverages > 1
    nFrames = size(data, 4)
    nBlocks = ceil(Int, nFrames / numAverages)

    rem(nFrames, numAverages) != 0 && (warn("numAverages no integer divisor of nFrames.
                     Last Block will be averaged over less than $numAverages Frames."))

    averagedData = zeros(eltype(data), size(data,1), size(data,2), size(data,3), nBlocks)

    p = Progress(nBlocks, 1, "Loading measurement from $(filepath(f)) ...")
    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nFrames) # ensure that modulo is taken into account

      tmp = data[:, :, :, index1:index2]

      averagedData[:, :, :, i] = mean(tmp,4)
      next!(p)
    end
    data = averagedData
  end

  return data
end

function getMeasurements(f::MPIFile; frames=1:acqNumFrames(f),
      loadas32bit=true, loadasreal=false, fourierTransform=true, kargs...)
  if experimentHasProcessing(f)
    data = getMeasurementsLowLevelFromProc(f; frames=frames, kargs...)
  else
    idx = measFGFrameIdx(f)
    data = getMeasurementsLowLevel(f; frames=idx[frames], kargs...)
  end

  data = loadas32bit ? map(Float32,data) : map(Float64,data)

  if fourierTransform
    data = rfft(data,1)
  end
  if loadasreal
    return returnasreal(data)
  else
    return data
  end
end

function getMeasurements(f::MPIFile, frequencies; loadasasreal=true, kargs...)
  data = getMeasurements(f, fourierTransform=true, loadasasreal=false; kargs...)
  data = reshape(data, size(data,1)*size(data,2), size(data,3), size(data,4))
  data = data[frequencies, :, :]
  if loadasreal
    return returnasreal(data)
  else
    return data
  end
end
