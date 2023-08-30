export filterFrequencies

"""
  filterFrequencies(f; kargs...) => Vector{Int64}

Supported keyword arguments:
* SNRThresh
* minFreq
* maxFreq
* recChannels
* sortBySNR
* numUsedFreqs
* stepsize
* maxMixingOrder
* sortByMixFactors
"""
function filterFrequencies(f::MPIFile; SNRThresh=-1, minFreq=0,
               maxFreq=rxBandwidth(f), recChannels=1:rxNumChannels(f),
               sortBySNR=false, numUsedFreqs=-1, stepsize=1, maxMixingOrder=-1,
               sortByMixFactors=false, numPeriodAverages=1, numPeriodGrouping=1)

  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  nReceivers = rxNumChannels(f)
  nPeriods = 1 #acqNumPeriodsPerFrame(f)
  freqIndices = collect(vec(CartesianIndices((nFreq, nReceivers))))

  minIdx = floor(Int, minFreq / rxBandwidth(f) * (nFreq-1) ) + 1
  maxIdx = ceil(Int, maxFreq / rxBandwidth(f) * (nFreq-1) ) + 1

  filterFrequenciesBySelection!(freqIndices, f)
  filterFrequenciesByChannel!(freqIndices, recChannels)
  filterFrequenciesByMinIdx!(freqIndices, minIdx)
  filterFrequenciesByMaxIdx!(freqIndices, maxIdx)

  #=
  if maxMixingOrder > 0
    if numPeriodGrouping == 1
      mf = mixingFactors(f)
      for l=1:size(mf,1)
        if mf[l,4] > maxMixingOrder || mf[l,4] > maxMixingOrder
          freqMask[(l-1)+1,recChannels] .= false
        end
      end
    else # This is a hack until we have a general solution

      numPatches = div(acqNumPeriodsPerFrame(f), numPeriodAverages)
      fBands = (collect(0:(rxNumFrequencies(f)-1)).-1) .* numPatches
      freqMaskMO = zeros(Bool,nFreq,nReceivers,nPeriods)

      for i=1:length(fBands)
        for y=-maxMixingOrder:maxMixingOrder
          index = fBands[i]+y
          if 1 <= index <= size(freqMask,1)
            freqMaskMO[index,recChannels] .= true
          end
        end
      end
      freqMask .&= freqMaskMO
    end
  end
  =#
  SNRAll = nothing

  if SNRThresh > 0 || numUsedFreqs > 0 || sortBySNR
    SNR = zeros(nFreq, nReceivers)
    idx = measIsFrequencySelection(f) ? measFrequencySelection(f) : idx = 1:nFreq

    SNRAll = calibSNR(f)
    if SNRAll != nothing
      SNR[idx,:] = SNRAll[:,:,1]
    end
  end

  if SNRThresh > 0 && numUsedFreqs > 0
    error("It is not possible to use SNRThresh and SNRFactorUsedFreq similtaneously")
  end

  if SNRThresh > 0 && SNRAll != nothing
    freqMask[ findall(vec(SNR) .< SNRThresh) ] .= false
  end

  if numUsedFreqs > 0 && SNRAll != nothing
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

  if sortBySNR && !sortByMixFactors && SNRAll != nothing
    SNR = vec(SNR[1:stepsize:nFreq,:,:])

    freq = freq[reverse(sortperm(SNR[freq]),dims=1)]
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

function filterFrequenciesBySelection!(indices::Vector{CartesianIndex}, f::MPIFile)
  if measIsFrequencySelection(f)
    return filterFrequenciesBySelection!(indices, measFrequencySelection(f))
  end
end
filterFrequenciesBySelection!(indices::Vector{CartesianIndex}, selection::Vector{Int64}) = filter!(x -> in(x[1], selection), indices)

filterFrequenciesByChannel!(indices::Vector{CartesianIndex}, channels) = filter!(x-> in(x[2], channels), indices)

function filterFrequenciesByMinFreq!(indices::Vector{CartesianIndex}, f::MPIFile, minFreq)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  minIdx = floor(Int, minFreq / rxBandwidth(f) * (nFreq-1) ) + 1
  return filterFrequenciesByMinIdx!(indices, minIdx)
end
filterFrequenciesByMinIdx!(indices::Vector{CartesianIndex}, minIdx) = minIdx > 0 ?  filter!(x -> x[1] > minFreq, freqIndices) : indices 


function filterFrequenciesByMaxFreq!(indices::Vector{CartesianIndex}, f::MPIFile, maxFreq)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  maxIdx = ceil(Int, maxFreq / rxBandwidth(f) * (nFreq-1) ) + 1
  return filterFrequenciesByMaxIdx!(indices, maxIdx)
end
filterFrequenciesByMaxIdx!(indices::Vector{CartesianIndex}, maxIdx) = filter!(x-> x[1] <= maxIdx, indices)

function sortFrequenciesBySNR!()
end



function sortFrequenciesByMixFactors()
end


function rowsToSubsampledRows(f::MPIFile, rows)
  if measIsFrequencySelection(f)
    # In this case we need to convert indices
    tmp = zeros(Int64, rxNumFrequencies(f), rxNumChannels(f) )
    idxAvailable = measFrequencySelection(f)
    for d=1:rxNumChannels(f)
      tmp[idxAvailable, d] = (1:length(idxAvailable)) .+ (d-1)*length(idxAvailable)
    end
    rows_ = vec(tmp)[rows]
    if findfirst(x -> x == 0, rows_) != nothing
      @error "Indices applied to systemMatrix are not available in the file"
    end
  else
    rows_ = rows
  end
  return rows_
end


