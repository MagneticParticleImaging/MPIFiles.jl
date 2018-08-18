export filterFrequencies

function filterFrequencies(f::MPIFile; SNRThresh=-1, minFreq=0,
               maxFreq=rxBandwidth(f), recChannels=1:rxNumChannels(f),
               sortBySNR=false, numUsedFreqs=-1, stepsize=1, maxMixingOrder=-1,
               sortByMixFactors=false)

  nFreq = rxNumFrequencies(f)
  nReceivers = rxNumChannels(f)
  nPeriods = 1 #acqNumPeriodsPerFrame(f)

  minIdx = round(Int, minFreq / rxBandwidth(f) * nFreq )
  maxIdx = round(Int, maxFreq / rxBandwidth(f) * nFreq )

  freqMask = zeros(Bool,nFreq,nReceivers,nPeriods)

  freqMask[:,recChannels,:] .= true

  if minIdx > 0
    freqMask[1:(minIdx-1),:,:] .= false
  end
  if maxIdx < nFreq
    freqMask[(maxIdx+1):end,:,:] .= false
  end

  if maxMixingOrder > 0
      mf = mixingFactors(f)
      for l=1:size(mf,1)
        if mf[l,4] > maxMixingOrder || mf[l,4] > maxMixingOrder
          freqMask[(l-1)+1,recChannels] = false
        end
      end
  end

  if SNRThresh > 0 || numUsedFreqs > 0 || sortBySNR
    SNR = calibSNR(f)[:,:,1]
  end

  if SNRThresh > 0 && numUsedFreqs > 0
    error("It is not possible to use SNRThresh and SNRFactorUsedFreq similtaneously")
  end

  if SNRThresh > 0
    freqMask[ findall(vec(SNR) .< SNRThresh) ] .= false
  end

  if numUsedFreqs > 0
    numFreqsAlreadyFalse = sum(!freqMask)
    numFreqsFalse = round(Int,length(freqMask)* (1-numUsedFreqs))
    S = sortperm(vec(SNR))

    l = 1
    j = 1
    while j<  (numFreqsFalse-numFreqsAlreadyFalse)
      if freqMask[S[l]] == true
        freqMask[S[l]] = false
        j += 1
      end
      l += 1
    end

  end


  if stepsize > 1
    freqStepsizeMask = zeros(Bool,nFreq, nReceivers, nPatches)
    freqStepsizeMask[1:stepsize:nFreq,:,:] = freqMask[1:stepsize:nFreq,:,:]
    freqMask = freqStepsizeMask
  end

  freq = findall( vec(freqMask) )

  if sortBySNR && !sortByMixFactors
    SNR = vec(SNR[1:stepsize:nFreq,:,:])

    freq = freq[flipdim(sortperm(SNR[freq]),1)]
  end

  if !sortBySNR && sortByMixFactors
    mfNorm = zeros(nFreq,nReceivers,nPeriods)
    mf = mixingFactors(f)
    for k=1:nFreq
      mfNorm[k,:,:] = norm(mf[k,1:3])
    end

    freq = freq[sortperm(mfNorm[freq])]
  end

  freq
end
