export filterFrequencies, sortFrequencies

function getCalibSNR(f::MPIFile; numPeriodGrouping = 1, stepsize = 1)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  nReceivers = rxNumChannels(f)
  SNR = zeros(nFreq, nReceivers)
  idx = measIsFrequencySelection(f) ? measFrequencySelection(f) : idx = 1:nFreq

  SNRAll = calibSNR(f)
  if SNRAll != nothing
    SNR[idx,:] = SNRAll[:,:,1]
  end

  SNR = SNR[1:stepsize:nFreq,:,:]
  return SNR
end

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
               maxFreq=rxBandwidth(f), recChannels=1:rxNumChannels(f),numUsedFreqs=-1, stepsize=1,
                maxMixingOrder=-1, numPeriodAverages=1, numPeriodGrouping=1, numSidebandFreqs = -1)

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

  if maxMixingOrder > 0
    if numPeriodGrouping == 1
      filterFrequenciesByMaxMixingOrder!(freqIndices, maxMixingOrder, f)
    else
      error("Can not apply max mixing order with a period grouping larger than one, found: $numPeriodGrouping")
    end
  end

  if numSidebandFreqs > 0 && numPeriodGrouping > 1
    filterFrequenciesByNumSidebandFreqs!(freqIndices, numSidebandFreqs, f, numPeriodGrouping = numPeriodGrouping)
  end

  SNRAll = nothing

  if SNRThresh > 0 || numUsedFreqs > 0
    SNR = zeros(nFreq, nReceivers)
    idx = measIsFrequencySelection(f) ? measFrequencySelection(f) : idx = 1:nFreq

    SNRAll = calibSNR(f)
    if SNRAll != nothing
      SNR[idx,:] = SNRAll[:,:,1]
    end
  end

  if SNRThresh > 0 && numUsedFreqs <= 0
    filterFrequenciesBySNRThresh!(freqIndices, SNRThresh, SNR)
  elseif  numUsedFreqs > 0 && SNRThresh <= 0
    filterFrequenciesByNumUsedFrequencies!(freqIndices, numUsedFreqs)
  elseif numUsedFreqs > 0 && SNRThresh > 0
    error("It is not possible to use SNRThresh and SNRFactorUsedFreq similtaneously")
  end


  if stepsize > 1
    filterFrequenciesByStepsize!(freqIndices, stepsize)
  end

  return freqIndices
end

export filterFrequenciesBySelection!
function filterFrequenciesBySelection!(indices::Vector{CartesianIndex{2}}, f::MPIFile)
  if measIsFrequencySelection(f)
    return filterFrequenciesBySelection!(indices, measFrequencySelection(f))
  end
end
filterFrequenciesBySelection!(indices::Vector{CartesianIndex{2}}, selection::Vector{Int64}) = filter!(x -> in(x[1], selection), indices)

export filterFrequenciesByChannel!
filterFrequenciesByChannel!(indices::Vector{CartesianIndex{2}}, channels) = filter!(x-> in(x[2], channels), indices)

export filterFrequenciesByMinFreq!
function filterFrequenciesByMinFreq!(indices::Vector{CartesianIndex{2}}, f::MPIFile, minFreq; numPeriodGrouping = 1)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  minIdx = floor(Int, minFreq / rxBandwidth(f) * (nFreq-1) ) + 1
  return filterFrequenciesByMinIdx!(indices, minIdx)
end
filterFrequenciesByMinIdx!(indices::Vector{CartesianIndex{2}}, minIdx) = minIdx > 0 ?  filter!(x -> x[1] > minIdx, indices) : indices 

export filterFrequenciesByMaxFreq!
function filterFrequenciesByMaxFreq!(indices::Vector{CartesianIndex{2}}, f::MPIFile, maxFreq; numPeriodGrouping = 1)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  maxIdx = ceil(Int, maxFreq / rxBandwidth(f) * (nFreq-1) ) + 1
  return filterFrequenciesByMaxIdx!(indices, maxIdx)
end
filterFrequenciesByMaxIdx!(indices::Vector{CartesianIndex{2}}, maxIdx) = filter!(x-> x[1] < maxIdx, indices)

export filterFrequenciesByMaxMixingOrder!
filterFrequenciesByMaxMixingOrder!(indices::Vector{CartesianIndex{2}}, maxMixingOrder, f::MPIFile) = filterFrequenciesByMaxMixingOrder!(indices, maxMixingOrder, mixingFactors(f))
filterFrequenciesByMaxMixingOrder!(indices::Vector{CartesianIndex{2}}, maxMixingOrder, mf::Matrix) = filter!(x-> mf[x[1], 4] <= maxMixingOrder, indices)

export filterFrequenciesByNumSidebandFreqs!
function filterFrequenciesByNumSidebandFreqs!(indices::Vector{CartesianIndex{2}}, numSidebandFreqs::Int64, f::MPIFile; numPeriodGrouping = 1)
  # https://en.wikipedia.org/wiki/Sideband
  # Because of period grouping we have more frequencies than in original data

  # Find "virtual" frequency indices of original frequencies
  fBands = (collect(0:(rxNumFrequencies(f))).-1) .* numPeriodGrouping

  delete = Int64[]
  for (i, cart) in enumerate(indices)
    freq = cart[1]
    # Check if there is no frequency that is at most numSideBandFreqs away from our original frequency
    if !any(fBand -> abs(fBand - freq) <= numSidebandFreqs, fBands)
      push!(delete, i)
    end
  end
  deleteat!(indices, delete)
end

export filterFrequenciesBySNRThresh!
function filterFrequenciesBySNRThresh!(indices::Vector{CartesianIndex{2}}, f::MPIFile, snrthresh; numPeriodGrouping = 1)
  SNR = getCalibSNR(f, numPeriodGrouping = numPeriodGrouping)
  return filterFrequenciesBySNRThresh!(indices, snrthresh, SNR)
end
filterFrequenciesBySNRThresh!(indices::Vector{CartesianIndex{2}}, SNRThresh, SNR) = filter!(x-> SNR[x] >= SNRThresh , indices)

export filterFrequenciesByNumUsedFrequencies!
function filterFrequenciesByNumUsedFrequencies!(indices::Vector{CartesianIndex{2}}, f::MPIFile, maxFreq)
  error("Not implemented")
end
filterFrequenciesByNumUsedFrequencies!(indices::Vector{CartesianIndex{2}}, maxIdx) = error("not implemented")
#=
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
=#

export filterFrequenciesByStepsize!
filterFrequenciesByStepsize!(indices::Vector{CartesianIndex{2}}, stepsize) = error("Not implemented")
#=
    freqStepsizeMask = zeros(Bool,nFreq, nReceivers, nPatches)
    freqStepsizeMask[1:stepsize:nFreq,:,:] = freqMask[1:stepsize:nFreq,:,:]
    freqMask = freqStepsizeMask
=#

function sortFrequencies!(indices::Vector{CartesianIndex{2}}, f::MPIFile; numPeriodGrouping = 1, stepsize = 1, sortBySNR = false, sortByMixFactors = false)
  if sortBySNR && !sortByMixFactors
    indices = sortFrequenciesBySNR(indices, f, numPeriodGrouping = numPeriodGrouping, stepsize = stepsize)
  elseif !sortBySNR && sortByMixFactors
    indices = sortFrequenciesByMixFactors(indices, f, numPeriodGrouping = numPeriodGrouping)
  else
    error("Can not apply multiple sorting algorithms to frequencies")
  end
  return indices
end

export sortFrequenciesBySNR
function sortFrequenciesBySNR(indices::Vector{CartesianIndex{2}}, f::MPIFile; numPeriodGrouping = 1, stepsize = 1)
  SNR = getCalibSNR(f, numPeriodGrouping = numPeriodGrouping, stepsize = stepsize)
  sortFrequenciesBySNR(indices, SNR)
end
sortFrequenciesBySNR(indices::Vector{CartesianIndex{2}}, SNR::AbstractArray) = indices[reverse(sortperm(SNR[indices]),dims=1)]

export sortFrequenciesByMixFactors
function sortFrequenciesByMixFactors(indices, f::MPIFile; numPeriodGrouping = 1)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  nReceivers = rxNumChannels(f)
  nPeriods = 1 #acqNumPeriodsPerFrame(f)
  mfNorm = zeros(nFreq,nReceivers,nPeriods)
  mf = mixingFactors(f)
  for k=1:nFreq
    mfNorm[k,:,:] = norm(mf[k,1:3])
  end
  sortFrequenciesByMixFactors(indices, mfNorm)
end
sortFrequenciesByMixFactors(indices::Vector{CartesianIndex{2}}, mfNorm::AbstractArray) = indices[sortperm(mfNorm[indices])]


function rowsToSubsampledRows(f::MPIFile, rows)
  if measIsFrequencySelection(f)
    # In this case we need to convert indices
    tmp = Array{Union{Nothing, CartesianIndex{2}}}(undef, rxNumFrequencies(f), rxNumChannels(f))
    idxAvailable = measFrequencySelection(f)
    for (i, idx) in enumerate(idxAvailable)
      for d=1:rxNumChannels(f)
        tmp[idx, d] = CartesianIndex(i, d)
      end
    end
    rows_ = tmp[rows]
    if any(isnothing, rows_)
      @error "Indices applied to systemMatrix are not available in the file"
    end
    return identity.(rows_) # Go from Vector{Union{Nothing, Index}} to Vector{Index}
  else
    rows_ = rows
  end
  return rows_
end


