####### Mixing factor functions #########

export mixFactorToFreq, mixFactorToFreqIdx, mixingFactors


"""
This function returns the multiplicator `k` corresponding to the frequency
given by the mixing factors `mx`, `my`, and `mz`. I.e. the integer value, for
which k*F = mx*fx + my*fy* + mz*fz, where F is the measurement cycle frequency
and fx, fy, and fz are the excitation frequencies for the x,y, and z channel.
"""
function mixFactorToFreq(b::MPIFile,mixFactors...)
  prefactors = calcPrefactors(b)
  if length(mixFactors) != length(prefactors)
    error("The MPIFile defines $(length(prefactors)) prefactors but you provided $(length(mixFactors)) mixing factors!")
  end
  k = sum(mixFactors.*prefactors)
  return k
end

"""
This function returns the index `freqidx` of the frequency given by the mixing
factors `mx`, `my`, and `mz` with respect to the frequency list
`freq = frequencies(bSF)`.
"""
mixFactorToFreqIdx(b::MPIFile,mx...) = mixFactorToFreq(b,mx...)+1


function calcPrefactors(dfStrength, divider, baseFreq, cycle)
  mask = collect((dfStrength[1,:,1] .>= 0.0000001))
  divider = vec(divider)
  #mxyz = round.(Int64,divider.*mask./gcd(divider.*mask))
  mxyz_ = baseFreq*cycle./divider
  mxyz = max.(1,round.(Int64,mxyz_.*mask))

  return Tuple(mxyz)
end

calcPrefactors(b::MPIFile) = calcPrefactors(dfStrength(b), dfDivider(b), dfBaseFrequency(b), dfCycle(b))

function calcPrefactors(f::MDFFileV2)
  if haskey(f.file, "/measurement/_manualPrefactors")
    prefactors = f.file["/measurement/_manualPrefactors"]
    return Tuple(prefactors[])
  else
    return calcPrefactors(dfStrength(f), dfDivider(f), dfBaseFrequency(f), dfCycle(f))
  end
end

"""
This function returns a lookup table with columns
`freqidx` | `mx` | `my` | `mz` | `mixingOrder`
for all frequencies in `freq = frequencies(bSF)`, where only the lowest order
mixing coefficients `mx`, `my`, and `mz` are listed.
"""
function mixingFactors(b::MPIFile; maxFactor=100)
  prefactors = calcPrefactors(b)
  return _mixingFactors(prefactors, rxNumFrequencies(b), maxFactor)
end

function _mixingFactors(prefactors::NTuple{Nt,Int}, numFreqs::Integer, N::Integer) where Nt
  mixingOrderList = zeros(Int, numFreqs, length(prefactors)+1)
  mixingOrderList[:,end] .= -1
  loop = CartesianIndices(ntuple(x->(-N:N), length(prefactors)))
  for mixFactors in loop
    mixFactorsTup = Tuple(mixFactors)
    k = sum(mixFactorsTup.*prefactors)
    mixOrder = mapreduce(abs,+,mixFactorsTup)
    if 1<=k<=numFreqs && (mixingOrderList[k,end]<0 || mixingOrderList[k,end]>=mixOrder)
          mixingOrderList[k,1:end-1] .= mixFactorsTup
          mixingOrderList[k,end] = mixOrder
        end
  end
  return mixingOrderList
end