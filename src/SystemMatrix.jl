export getSystemMatrix, getSystemMatrixReshaped, calculateSystemMatrixSNR

function converttoreal(S::AbstractArray{Complex{T},2}) where {T}
  N = size(S,1)
  M = size(S,2)
  S = reshape(reinterpret(T,vec(S)),(2*N,M))
  p = Progress(M, 1, "Converting system matrix to real...")
  for l=1:M
    tmp = S[:,l]
    S[1:N,l] = tmp[1:2:end]
    S[N+1:end,l] = tmp[2:2:end]
    next!(p)
  end
  return reshape(S,(N,2*M))
end

function getSystemMatrix(f::MPIFile,
           frequencies=1:rxNumFrequencies(f)*rxNumChannels(f)*acqNumPeriodsPerFrame(f);
                         bgCorrection=false, loadasreal=false,
                         kargs...)
  #if measIsTransposed(f) && measIsFourierTransformed(f)
  data = systemMatrix(f, frequencies, bgCorrection)
  #else
  #  error("TODO: implement making a SF using getMeasurement")
  #end
  S = map(ComplexF32, data)

  if loadasreal
    return converttoreal(S)
  else
    return S
  end
end

function getSystemMatrixReshaped(f::MPIFile; kargs...)
  return reshape(getSystemMatrix(f;kargs...),:,rxNumFrequencies(f),
                            rxNumChannels(f),acqNumPeriodsPerFrame(f))
end

#function calculateSystemMatrixSNR(f::MPIFile)
#
#end

function calculateSystemMatrixSNR(f::MPIFile, S::Array)
  SNR = zeros(rxNumFrequencies(f),rxNumChannels(f),acqNumPeriodsPerFrame(f))
  for j=1:acqNumPeriodsPerFrame(f)
    for r=1:rxNumChannels(f)
      for k=1:rxNumFrequencies(f)
        meanBG = mean(S[(acqNumFGFrames(f)+1):end,k,r,j])
        signal = maximum(abs.(S[1:acqNumFGFrames(f),k,r,j].-meanBG))
        noise = mean(abs.(S[(acqNumFGFrames(f)+1):end,k,r,j].-meanBG))
        SNR[k,r,j] = signal / noise
      end
    end
  end
  return SNR
end
