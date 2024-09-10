export getMeasurements, getMeasurementsFD

function measDataConv(f::MPIFile, args...)
  data = measDataTD(f, args...)

  if eltype(data) <: Integer
    data = map(Float32, data)
  end
  a = rxDataConversionFactor(f)
  if !isnothing(a)
    for d=1:size(data,2)
      slice = view(data,:,d,:,:)
      rmul!(slice, a[1,d])
      slice .+= a[2,d]
    end
  end
  return data
end

hannWindow(M) = (1.0 .- cos.(2*π/(M-1)*(0:(M-1))))/(M-1)*M

function measDataSpectralLeakageCorrectedSinglePatch(f::MPIFile, frames)
  @debug "Apply Spectral Cleaning"
  numTimePoints = rxNumSamplingPoints(f)
  numReceivers = rxNumChannels(f)
  numFrames = acqNumFrames(f)

  data = zeros(Float32, numTimePoints, numReceivers, 1, length(frames))

  window3 = hannWindow(numTimePoints*3)
  window2 = hannWindow(numTimePoints*2)

  for (i,fr) in enumerate(frames)
    for r in 1:numReceivers
      if fr==1
        tmp = measDataConv(f, fr:fr+1, 1, r)
        data[:,r,1,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                          +  tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      elseif fr==numFrames
        tmp = measDataConv(f, fr-1:fr, 1, r)
        data[:,r,1,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                          +    tmp[:,1,1,2] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      else
        tmp = measDataConv(f, fr-1:fr+1, 1, r)
        data[:,r,1,i] = 1/3 * (tmp[:,1,1,1] .* window3[1:numTimePoints]
                          +  tmp[:,1,1,2] .* window3[1+numTimePoints:2*numTimePoints]
                          +  tmp[:,1,1,3] .* window3[1+2*numTimePoints:3*numTimePoints]
                          );
      end
    end
  end
  return data
end

function measDataSpectralLeakageCorrectedMultiPatch(f::MPIFile, frames, periods)
  @debug "Apply Spectral Cleaning"
  numTimePoints = rxNumSamplingPoints(f)
  numReceivers = rxNumChannels(f)
  numFrames = acqNumFrames(f)
  numPeriods = acqNumPeriodsPerFrame(f)

  data = zeros(Float32, numTimePoints, numReceivers, length(periods), length(frames))

  window3 = hannWindow(numTimePoints*3)
  window2 = hannWindow(numTimePoints*2)

  for (i,fr) in enumerate(frames)
    for (p,pe) in enumerate(periods)
      for r in 1:numReceivers
        if pe == 1
          tmp = measDataConv(f, fr, pe:pe+1, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +  tmp[:,1,2,1] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        elseif pe == numPeriods
          tmp = measDataConv(f, fr, pe-1:pe, r)
          data[:,r,p,i] = 1/2 * (tmp[:,1,1,1] .* window2[1:numTimePoints]
                            +    tmp[:,1,2,1] .* window2[1+numTimePoints:2*numTimePoints]
                            );
        else
          tmp = measDataConv(f, fr, pe-1:pe+1, r)
          data[:,r,p,i] = 1/3 * (tmp[:,1,1,1] .* window3[1:numTimePoints]
                            +  tmp[:,1,2,1] .* window3[1+numTimePoints:2*numTimePoints]
                            +  tmp[:,1,3,1] .* window3[1+2*numTimePoints:3*numTimePoints]
                            );
        end
      end
    end
  end
  return data
end

function measDataSpectralLeakageCorrected(f::MPIFile, frames, periods)
  if acqNumPeriodsPerFrame(f) == 1
    return measDataSpectralLeakageCorrectedSinglePatch(f, frames)
  else
    return measDataSpectralLeakageCorrectedMultiPatch(f, frames, periods)
  end
end

function measDataLowLevel(f::MPIFile, args...; spectralLeakageCorrection=false)
  if !spectralLeakageCorrection || measIsSpectralLeakageCorrected(f) ||
      acqNumFrames(f) == 1
     tmp = measDataConv(f, args...)
  else
     tmp = measDataSpectralLeakageCorrected(f, args...)
  end
  return tmp
end

function returnasreal(u::AbstractArray{Complex{T}}) where {T}
  return copy(reshape(reinterpret(T,vec(u)),tuple(size(u,1)*2,size(u)[2:end]...)))
end
returnasreal(u::AbstractArray{T}) where {T<:Real} = u

function getAveragedMeasurements(f::MPIFile; frames=1:acqNumFrames(f),
            numAverages=1,  periods=1:acqNumPeriodsPerFrame(f),
            averagePeriodsPerPatch=false, numPeriodAverages=1, fixDistortions=false, kargs...)

  @debug "frequency and frame selection" rxNumSamplingPoints(f) rxNumChannels(f) acqNumFrames(f)

  if numAverages == 1
    data = measDataLowLevel(f, frames, periods; kargs...)
  else
    nFrames = length(frames)
    nBlocks = ceil(Int, nFrames / numAverages)

    if rem(nFrames, numAverages) != 0
      @warn "numAverages no integer divisor of nFrames ($nFrames).
              Last Block will be averaged over less than $numAverages Frames."
    end

    data = zeros(Float32, rxNumSamplingPoints(f), rxNumChannels(f), acqNumPeriodsPerFrame(f), nBlocks)

    for i = 1:nBlocks
      index1 = 1 + (i-1)*numAverages
      index2 = min( index1 + numAverages-1, nFrames) # ensure that modulo is taken into account

      tmp = measDataLowLevel(f, frames[index1:index2], periods; kargs...)
      data[:,:,:,i] = mean(tmp,dims=4)
    end
  end

  if fixDistortions
    detectAndFixDistortions!(data, 0.3)
  end

  if numPeriodAverages > 1 && averagePeriodsPerPatch
    error("getAveragedMeasurements: not possible to combine numPeriodAverages and averagePeriodsPerPatch")
  end

  if averagePeriodsPerPatch
    if periods != 1:acqNumPeriodsPerFrame(f)
      error("Option averagePeriodsPerPatch can only be used when all periods are selected")
    end
    data_ = reshape(data, rxNumSamplingPoints(f), rxNumChannels(f),
                          acqNumPeriodsPerPatch(f), acqNumPatches(f), size(data,4))
    dataAv = mean(data_, dims=3)

    return reshape(dataAv, rxNumSamplingPoints(f), rxNumChannels(f), acqNumPatches(f), size(data,4))
  elseif numPeriodAverages > 1
    newNumPeriods = div(acqNumPeriodsPerFrame(f), numPeriodAverages)
    data_ = reshape(data, rxNumSamplingPoints(f), rxNumChannels(f), numPeriodAverages, newNumPeriods, size(data,4))
    dataAv = mean(data_, dims=3)
    return reshape(dataAv, rxNumSamplingPoints(f), rxNumChannels(f), newNumPeriods, size(data,4))    
  else
    return data
  end
end


"""
  getMeasurements(f, [neglectBGFrames]; kargs...) => Array{Float32,4}

Load the measurement data in time domain

Supported keyword arguments:
* frames
* bgCorrection
* interpolateBG
* tfCorrection
* sortFrames
* numAverages
* numPeriodAverages
* spectralLeakageCorrection
"""
function getMeasurements(f::MPIFile, neglectBGFrames=true;
      frames=neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f)),
      bgCorrection=false, bgFrames = 1:length(measBGFrameIdx(f)), interpolateBG=false, tfCorrection=rxHasTransferFunction(f),
      sortFrames=false, numAverages=1, numPeriodGrouping=1, kargs...)

  if neglectBGFrames
    idx = measFGFrameIdx(f)

    data = getAveragedMeasurements(f; frames=idx[frames],
                                      numAverages=numAverages, kargs...)
    
    idxBG = measBGFrameIdx(f)[bgFrames]
    hasBGFrames = length(idxBG) > 0
    if bgCorrection && !hasBGFrames
      @warn "Background correction was selected but there are no background frames in the file."
    elseif bgCorrection && hasBGFrames
      @debug "Applying bg correction ..."
      
      dataBG = getAveragedMeasurements(f; frames=idxBG, kargs...)
      if interpolateBG
        blockLen = measBGFrameBlockLengths(measIsBGFrame(f))
        st = 1
        for j=1:length(blockLen)
          dataBG[:,:,:,st:st+blockLen[j]-1] .=
               mean(dataBG[:,:,:,st:st+blockLen[j]-1], dims=4)
          st += blockLen[j]
        end

        dataBGInterp = interpolate(dataBG,
          (NoInterp(), NoInterp(), NoInterp(), BSpline(Linear()))) #OnCell?

        origIndex = idx[frames]
        M = size(data,4)
        K = size(dataBG,4)
        N = M + K
        for m=1:M
          alpha = (origIndex[m]-1)/(N-1)*(K-1)+1
          for k1=1:size(data,1)
            for k2=1:size(data,2)
              for k3=1:size(data,3)
                data[k1,k2,k3,m] -= dataBGInterp(k1,k2,k3,alpha)
              end
            end
          end
        end
      else
        data[:,:,:,:] .-= mean(dataBG, dims=4)
      end
    end

    if sortFrames
      if calibIsMeanderingGrid(f)
        data[:,:,:,:] = data[:,:,:,meanderingFramePermutation(f)]
      end
    end
  else
    if sortFrames
      permJoint = fullFramePermutation(f)
      data = getAveragedMeasurements(f; frames=permJoint, numAverages=numAverages, kargs...)
    else
      data = getAveragedMeasurements(f; frames=frames, numAverages=numAverages, kargs...)
    end

    if bgCorrection && hasBGFrames
      dataBG = getAveragedMeasurements(f; frames=idxBG, kargs...)

      data[:,:,:,:] .-= mean(dataBG, dims=4)
    end
  end

  if numPeriodGrouping > 1
    tmp = permutedims(data, (1,3,2,4))

    if mod(size(tmp,2),numPeriodGrouping) != 0
      error("Periods cannot be grouped because $(size(tmp,2)) cannot be divided by $numPeriodGrouping")
    end

    tmp2 = reshape(tmp, size(tmp,1)*numPeriodGrouping, div(size(tmp,2),numPeriodGrouping),
                        size(tmp,3), size(tmp,4) )
    data = permutedims(tmp2, (1,3,2,4))
  end

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    inductionFactor = rxInductionFactor(f)
    if isnothing(tf)
      error("No transfer function available in file, please use tfCorrection=false")
    end
    if isnothing(inductionFactor)
      @warn "The file is missing the induction factor. The induction factor will be set to 1."
      inductionFactor = ones(Float64, rxNumChannels(f))
    end

    J = size(data,1)
    dataF = rfft(data, 1)
    dataF ./= tf
    dataF[isnan.(dataF)] .= zero(eltype(dataF))

    @warn "This measurement has been corrected with a Transfer Function. Name of TF: $(rxTransferFunctionFileName(f))"
    if !isnothing(inductionFactor)
        K = minimum([length(inductionFactor), size(dataF, 2)])
        if K != length(inductionFactor)
          @warn "The amount of channels in the data and the TF does not match. Using lowest value. Please check closely if the TF is applied to wrong channels."
        end

       	for k=1:K
       		dataF[:,k,:,:] ./= inductionFactor[k]
       	end
    end
    data = irfft(dataF,J,1)
  end

  return map(Float32,data)
end


"""
  getMeasurementsFD(f, [neglectBGFrames]; kargs...) => Array{ComplexF32,4}

Load the measurement data in frequency domain

Supported keyword arguments:
* frames
* bgCorrection
* interpolateBG
* tfCorrection
* sortFrames
* numAverages
* spectralLeakageCorrection
* loadasreal
* transposed
* frequencies
"""
function getMeasurementsFD(f::MPIFile, args...;
      loadasreal=false, transposed=false, frequencies=nothing,
      tfCorrection=rxHasTransferFunction(f), kargs...)

  data = getMeasurements(f, args..., tfCorrection=false; kargs...)
  
  data = rfft(data, 1)

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    inductionFactor = rxInductionFactor(f)
    if isnothing(tf)
      error("No transfer function available in file, please use tfCorrection=false")
    end
    if isnothing(inductionFactor)
      @warn "The file is missing the induction factor. The induction factor will be set to 1."
      inductionFactor = ones(Float64, rxNumChannels(f))
    end

    data ./= tf

    @warn "This measurement has been corrected with a Transfer Function. Name of TF: $(rxTransferFunctionFileName(f))"
    if !isnothing(inductionFactor)
        K = minimum([length(inductionFactor), size(data, 2)])
        if K != length(inductionFactor)
          @warn "The amount of channels in the data and the TF does not match. Using lowest value. Please check closely if the TF is applied to wrong channels."
        end

       	for k=1:K
       		data[:,k,:,:] ./= inductionFactor[k]
       	end
    end
  end

  if !isnothing(frequencies)
    # here we merge frequencies and channels
    data = data[frequencies, :, :]
  end

  if transposed
    if !isnothing(frequencies)
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

function spectralLeakageCorrectedData(dataIn)
  @debug "Apply Spectral Cleaning"

  numTimePoints = size(dataIn,1)
  numFrames = size(dataIn,2)

  dataOut = zeros(Float32, numTimePoints, numFrames)

  window3 = hannWindow(numTimePoints*3)
  window2 = hannWindow(numTimePoints*2)

  for (i,fr) in enumerate(collect(1:numFrames))
      if fr==1
        dataOut[:,i] = 1/2 * (dataIn[:,fr] .* window2[1:numTimePoints]
                          +  dataIn[:,fr+1] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      elseif fr==numFrames
        dataOut[:,i] = 1/2 * (dataIn[:,fr-1] .* window2[1:numTimePoints]
                          +    dataIn[:,fr] .* window2[1+numTimePoints:2*numTimePoints]
                          );
      else
        dataOut[:,i] = 1/3 * (dataIn[:,fr-1] .* window3[1:numTimePoints]
                          +  dataIn[:,fr] .* window3[1+numTimePoints:2*numTimePoints]
                          +  dataIn[:,fr+1] .* window3[1+2*numTimePoints:3*numTimePoints]
                          );
      end
    end
  return dataOut
end


function detectAndFixDistortions!(data::AbstractArray{T,4}, thresh) where T
  for fr=1:size(data, 4)
    peaks = Int[]
    for p = 1:size(data, 3)
      if maximum(abs.(data[:,1,p,fr])) > thresh
        push!(peaks, p)
      end
    end
    if !isempty(peaks)
     @show peaks
    end
    
    for q in peaks
      qbegin = q - 1
      qend = q + 1
      while qbegin > 0 && qbegin ∈ peaks
        qbegin -= 1
      end
      while qend <= size(data,3) && qend ∈ peaks
        qend += 1
      end
      if qbegin >=1 && qend <= size(data,3)
        α = (q - qbegin)/(qend - qbegin)
        data[:,:,q,fr] = α * data[:,:,qbegin,fr] + (1 - α) * data[:,:,qend,fr]
      end
    end
  end
  return
end