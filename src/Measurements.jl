export getMeasurements, getMeasurementsLowLevel


function measDataTDPeriodsConv(f::MPIFile, args...)
  data = measDataTDPeriods(f, args...)

  data_ = eltype(data) <: Integer ? map(Float32, data) : data

  a = rxDataConversionFactor(f)
  if a!=nothing
    for d=1:size(data_,2)
      slice = view(data_,:,d,:)
      scale!(slice, a[1,d])
      slice .+= a[2,d]
    end
  end
  return data_
end

hannWindow(m,M) = (1.-cos.(2*π/(M-1)*(m-1)))/(M-1)*M # TK: not sure about (M-1)*M ...

function measDataTDPeriodsSLCorr(f::MPIFile, periods, receivers)
  numTimePoints = rxNumSamplingPoints(f)
  numReceivers = rxNumChannels(f)
  numPeriods = acqNumPeriods(f)

  data = zeros(Float32, numTimePoints, length(receivers), length(periods))
  window3 = hannWindow(1:(numTimePoints*3),numTimePoints*3)
  window2 = hannWindow(1:(numTimePoints*2),numTimePoints*2)

  for (p,pe) in enumerate(periods)
    for (r,re) in enumerate(receivers)
      if pe == 1
        tmp = measDataTDPeriodsConv(f, pe:pe+1, re)
        data[:,r,p,i] = 1/2 * (tmp[:,1,1] .* window2[1:numTimePoints]
                          +  tmp[:,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      elseif pe == numPeriods
        tmp = measDataTDPeriodsConv(f, pe-1:pe, re)
        data[:,r,p,i] = 1/2 * (tmp[:,1,1] .* window2[1:numTimePoints]
                          +    tmp[:,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      else
        tmp = measDataTDPeriodsConv(f, pe-1:pe+1, re)
        data[:,r,p,i] = 1/3 * (tmp[:,1,1] .* window3[1:numTimePoints]
                          +  tmp[:,1,2] .* window3[1+numTimePoints:2*numTimePoints]
                          +  tmp[:,1,3] .* window3[1+2*numTimePoints:3*numTimePoints]
                          );
      end
    end
  end
  return data
end

function measDataLowLevel(f::MPIFile, args...; spectralLeakageCorrection=false )
  if !spectralLeakageCorrection || measIsSpectralLeakageCorrected(f) ||
    acqNumPeriods(f) == 1
    tmp = measDataTDPeriodsConv(f, args...)
  else
    tmp = measDataTDPeriodsSLCorr(f, args...)
  end
end


function returnasreal{T}(u::AbstractArray{Complex{T}})
  return reinterpret(T,u,tuple(size(u,1)*2,size(u)[2:end]...))
end
returnasreal{T<:Real}(u::AbstractArray{T}) = u

function getAveragedMeasurementsPeriods(f::MPIFile; periods=1:acqNumPeriods(f),
            numAverages=1, verbose = false, receivers=1:rxNumChannels(f),
            kargs...)

  if numAverages == 1
    data = measDataLowLevel(f, periods, receivers; kargs...)
  else
    nPeriods = length(periods)
    nBlocks = ceil(Int, nPeriods / numAverages)

    rem(nPeriods, numAverages) != 0 && (warn("numAverages no integer divisor of nPeriods.
              Last Block will be averaged over less than $numAverages Frames."))

    data = zeros(Float32, rxNumSamplingPoints(f), length(receivers), nBlocks)

    p = Progress(nBlocks, 1, "Loading measurement from $(filepath(f)) ...")
    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nPeriods) # ensure that modulo is taken into account

      tmp = measDataLowLevel(f, periods[index1:index2], receivers; kargs...)

      data[:,:,i] = mean(tmp,3)
      next!(p)
    end
  end
  return data
end

#function getAveragedMeasurements(f::MPIFile; frames=1:acqNumFrames(f), kargs...)
#  periods = ...
#  data = getAveragedMeasurementsPeriods(f, periods=periods; kargs...)
#  return reshape(data, size(data,1), size(data,2), :, length(frames))
#end


function getAveragedMeasurements(f::MPIFile; frames=1:acqNumFrames(f),
            numAverages=1,  verbose = false, periods=1:acqNumPeriodsPerFrame(f),
            spectralLeakageCorrection=false)

  verbose && println( rxNumSamplingPoints(f), " ",
                      rxNumChannels(f), " ", acqNumFrames(f), )

  if numAverages == 1
    data = measDataLowLevel(f, periods, spectralLeakageCorrection=spectralLeakageCorrection)
    data = reshape(data, size(data,1), size(data,2), :, length(frames))
  else
    nFrames = length(frames)
    nBlocks = ceil(Int, nFrames / numAverages)

    rem(nFrames, numAverages) != 0 && (warn("numAverages no integer divisor of nFrames.
              Last Block will be averaged over less than $numAverages Frames."))

    data = zeros(Float32, rxNumSamplingPoints(f), rxNumChannels(f), acqNumPeriodsPerFrame(f), nBlocks)

    p = Progress(nBlocks, 1, "Loading measurement from $(filepath(f)) ...")
    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nFrames) # ensure that modulo is taken into account

      tmp = measDataLowLevel(f, frames[index1:index2], periods, spectralLeakageCorrection=spectralLeakageCorrection)
      data[:,:,:,i] = mean(tmp,4)
      next!(p)
    end
  end
  return data
end


function getMeasurements(f::MPIFile, neglectBGFrames=true; periods=1:acqNumPeriods(f),
      loadasreal=false, fourierTransform=false, receivers=1:rxNumChannels(f),
      transposed=false, bgCorrection=false, frequencies=nothing,
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

  if fourierTransform || (frequencies != nothing)
    data = rfft(data, 1)
  end

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    if fourierTransform || (frequencies != nothing)
      data ./= tf
    else
      dim = 1
      J = size(data,dim)
      dataF = rfft(data, dim)
      dataF ./= tf
      data = irfft(dataF,J,dim)
    end
  end

  if frequencies != nothing
    # here we merge frequencies and channels
    data = reshape(data, size(data,1)*size(data,2), size(data,3), size(data,4))
    data = data[frequencies, :, :]
  end

  if transposed
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
