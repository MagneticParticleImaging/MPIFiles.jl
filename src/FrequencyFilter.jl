export filterFrequencies, sortFrequencies

function getCalibSNR(f::MPIFile; numPeriodGrouping = 1, stepsize = 1)
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  nReceivers = rxNumChannels(f)
  SNR = zeros(nFreq, nReceivers)
  idx = measIsFrequencySelection(f) ? measFrequencySelection(f) : idx = 1:nFreq

  SNRAll = calibSNR(f)
  if !isnothing(SNRAll)
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
* numPeriodGrouping
* numSideBandFreqs
* stopBands
"""
function filterFrequencies(f::MPIFile; SNRThresh=-1, minFreq=0,
               maxFreq=rxBandwidth(f), recChannels=1:rxNumChannels(f),numUsedFreqs=-1, stepsize=1,
                maxMixingOrder=-1, numPeriodAverages=1, numPeriodGrouping=1, numSidebandFreqs = -1, stopBands = nothing)

  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  nReceivers = rxNumChannels(f)
  nPeriods = 1 #acqNumPeriodsPerFrame(f)
  freqs = measIsFrequencySelection(f) ? measFrequencySelection(f) : 1:nFreq

  if numPeriodGrouping == 1
    freqIndices = vec([CartesianIndex{2}(i, j) for i in freqs, j in intersect(1:nReceivers, recChannels)])
  else
    freqIndices = collect(vec(CartesianIndices((nFreq, nReceivers))))
    filterFrequenciesBySelection!(freqIndices, f) 
    filterFrequenciesByChannel!(freqIndices, recChannels)
  end  
  if minFreq > 0
    filterFrequenciesByMinFreq!(freqIndices, f, minFreq)
  end
  
  if maxFreq < rxBandwidth(f)
    filterFrequenciesByMaxFreq!(freqIndices, f, maxFreq)
  end

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
    SNR = zeros(Float64, nFreq, nReceivers)

    SNRAll = calibSNR(f)
    if !isnothing(SNRAll)
      SNR[freqs, :] = SNRAll[:, :, 1]
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

  if !isnothing(stopBands)
    filterFrequenciesByStopBands!(freqIndices, f, stopBands; numPeriodGrouping = numPeriodGrouping)
  end

  return collect(vec(freqIndices))
end

export filterFrequenciesBySelection!
function filterFrequenciesBySelection!(indices, f::MPIFile)
  if measIsFrequencySelection(f)
    return filterFrequenciesBySelection!(indices, measFrequencySelection(f))
  end
end
filterFrequenciesBySelection!(indices, selection::Vector{Int64}, sorted = issorted(selection)) = filter!(x -> sorted ? insorted(x[1], selection) : in(x[1], selection), indices)

export filterFrequenciesByChannel!
filterFrequenciesByChannel!(indices, channels, sorted = issorted(channels)) = filter!(x-> sorted ? insorted(x[2], channels) : in(x[2], channels), indices)

export filterFrequenciesByMinFreq!
function filterFrequenciesByMinFreq!(indices, f::MPIFile, minFreq; numPeriodGrouping = 1)
  minIdx = searchsortedfirst(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), minFreq)
  return filterFrequenciesByMinIdx!(indices, minIdx)
end
filterFrequenciesByMinIdx!(indices, minIdx) = minIdx > 0 ?  filter!(x -> x[1] >= minIdx, indices) : indices 

export filterFrequenciesByMaxFreq!
function filterFrequenciesByMaxFreq!(indices, f::MPIFile, maxFreq; numPeriodGrouping = 1)
  maxIdx = searchsortedlast(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), maxFreq)
  return filterFrequenciesByMaxIdx!(indices, maxIdx)
end
filterFrequenciesByMaxIdx!(indices, maxIdx) = filter!(x-> x[1] <= maxIdx, indices)

export filterFrequenciesByMaxMixingOrder!
filterFrequenciesByMaxMixingOrder!(indices, maxMixingOrder, f::MPIFile) = filterFrequenciesByMaxMixingOrder!(indices, maxMixingOrder, mixingFactors(f))
filterFrequenciesByMaxMixingOrder!(indices, maxMixingOrder, mf::Matrix) = filter!(x-> 0 <= mf[x[1], 4] <= maxMixingOrder, indices)

export filterFrequenciesByNumSidebandFreqs!
function filterFrequenciesByNumSidebandFreqs!(indices, numSidebandFreqs::Int64, f::MPIFile; numPeriodGrouping = 1)
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
function filterFrequenciesBySNRThresh!(indices, f::MPIFile, snrthresh::T; numPeriodGrouping = 1) where T <: Real
  SNR = getCalibSNR(f, numPeriodGrouping = numPeriodGrouping)
  return filterFrequenciesBySNRThresh!(indices, snrthresh, SNR)
end
filterFrequenciesBySNRThresh!(indices, SNRThresh::T, SNR::Matrix) where T <: Real = filter!(x-> SNR[x] >= SNRThresh , indices)
filterFrequenciesBySNRThresh!(indices, SNRThresh::T, SNR::Dict{CartesianIndex{2}, Float64}) where T <: Real = filter!(x-> get(SNR, x, 0.0) >= SNRThresh , indices)


export filterFrequenciesByNumUsedFrequencies!
function filterFrequenciesByNumUsedFrequencies!(indices, f::MPIFile, maxFreq)
  error("Not implemented")
end
filterFrequenciesByNumUsedFrequencies!(indices, maxIdx) = error("not implemented")
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
function filterFrequenciesByStepsize!(indices, stepsize)
  stepIndices = 1:stepsize:maximum(map(x -> x[1], indices))
  filter!(x -> insorted(x[1], stepIndices), indices)
end

export filterFrequenciesByStopBands!, filterFrequenciesByStopBand!
function filterFrequenciesByStopBands!(indices, f::MPIFile, stopBands::Vector; kwargs...)
  for stopBand in stopBands
    filterFrequenciesByStopBand!(indices, f, stopBand; kwargs...)
  end
end
filterFrequenciesByStopBands!(indices, f::MPIFile, stopBand::Union{Vector{Int64}, UnitRange, NTuple}; kwargs...) = filterFrequenciesByStopBands!(indices, f, [stopBand]; kwargs...)
function filterFrequenciesByStopBand!(indices, f::MPIFile, stopBand::Vector{Int64}; numPeriodGrouping = 1)
  if length(stopBand) != 2
    error("Stop band are only defined for a start and stop value. Found $(length(stopBand)) values")
  end
  nFreq = rxNumFrequencies(f, numPeriodGrouping)
  minIdx = searchsortedlast(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), first(stopBand))
  maxIdx = searchsortedlast(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), last(stopBand))
  return filterFrequenciesByStopBand!(indices, minIdx, maxIdx)
end
function filterFrequenciesByStopBand!(indices, f::MPIFile, stopBand::Union{UnitRange, NTuple{2, Int64}}; numPeriodGrouping = 1)
  minIdx = searchsortedlast(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), first(stopBand))
  maxIdx = searchsortedlast(rfftfreq(rxNumSamplingPoints(f)*numPeriodGrouping, 2rxBandwidth(f)), last(stopBand))
  return filterFrequenciesByStopBand!(indices, minIdx, maxIdx)
end
filterFrequenciesByStopBand!(indices, minIdx, maxIdx) = filter!(x-> x[1] < minIdx || x[1] > maxIdx, indices)

function sortFrequencies(indices, f::MPIFile; numPeriodGrouping = 1, stepsize = 1, sortBySNR = false, sortByMixFactors = false)
  if sortBySNR && !sortByMixFactors
    indices = sortFrequenciesBySNR(indices, f, numPeriodGrouping = numPeriodGrouping, stepsize = stepsize)
  elseif !sortBySNR && sortByMixFactors
    indices = sortFrequenciesByMixFactors(indices, f, numPeriodGrouping = numPeriodGrouping)
  elseif sortBySNR && sortByMixFactors
    error("Can not apply multiple sorting algorithms to frequencies")
  end
  return indices
end

export sortFrequenciesBySNR
function sortFrequenciesBySNR(indices, f::MPIFile; numPeriodGrouping = 1, stepsize = 1)
  SNR = getCalibSNR(f, numPeriodGrouping = numPeriodGrouping, stepsize = stepsize)
  sortFrequenciesBySNR(indices, SNR)
end
sortFrequenciesBySNR(indices, SNR::AbstractArray) = indices[reverse(sortperm(SNR[indices]),dims=1)]

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
sortFrequenciesByMixFactors(indices, mfNorm::AbstractArray) = indices[sortperm(mfNorm[indices])]


function rowsToSubsampledRows(f::MPIFile, rows)
  if measIsFrequencySelection(f)
    # In this case we need to convert indices
    tmp = Array{Union{Nothing, CartesianIndex{2}}}(undef, rxNumFrequencies(f), rxNumChannels(f))
    idxAvailable = measFrequencySelection(f)
    for (i, idx) in enumerate(idxAvailable)
      for d=1:size(tmp, 2)
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


