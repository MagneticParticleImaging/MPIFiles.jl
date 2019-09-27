export getSystemMatrix, getSystemMatrixReshaped, calculateSystemMatrixSNR

"""
  getSystemMatrix(f, [neglectBGFrames]; kargs...) => Array{ComplexF32,4}

Load the system matrix in frequency domain

Supported keyword arguments:
* frequencies
* bgCorrection
* loadasreal
"""
function getSystemMatrix(f::MPIFile,
           frequencies=1:rxNumFrequencies(f)*rxNumChannels(f);
                         bgCorrection=false, loadasreal=false,
                         tfCorrection=rxHasTransferFunction(f), kargs...)

  data = systemMatrix(f, frequencies, bgCorrection)

  S = map(ComplexF32, data)

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    if tf != nothing
      _corrTFSF(S,tf[rowsToSubsampledRows(f,frequencies)])
    else
      error("TF not available")
    end
  end

  if loadasreal
    return converttoreal(S)
  else
    return S
  end
end

function _corrTFSF(S,tf)
  for l=1:size(S,2)
    S[:,l] ./= tf[l]
  end
end

function getSystemMatrixReshaped(f::MPIFile; kargs...)
  return reshape(getSystemMatrix(f;kargs...),:,rxNumFrequencies(f),
                            rxNumChannels(f),acqNumPeriodsPerFrame(f))
end

function calculateSystemMatrixSNR(f::MPIFile)
  data = systemMatrixWithBG(f)
  SNR = calculateSystemMatrixSNR(f, data)
  return SNR
end

function calculateSystemMatrixSNR(f::MPIFile, S::Array)
  J = acqNumPeriodsPerFrame(f)
  R = rxNumChannels(f)
  K = rxNumFrequencies(f)
  N = acqNumFGFrames(f)

  SNR = zeros(K, R, J)

  calculateSystemMatrixSNRInner(S, SNR, J, R, K, N)
  return SNR
end

function calculateSystemMatrixSNRInner(S, SNR, J, R, K, N)
  for j=1:J
    for r=1:R
      for k=1:K
        SBG = S[(N+1):end,k,r,j]
        SFG = S[1:N,k,r,j]
        diffBG = diff(SBG)
        meanBG = mean(SBG)
        signal = maximum(abs.(SFG.-meanBG))
        #noise = mean(abs.(SFG.-meanBG))
        noise = mean(abs.(diffBG))
        SNR[k,r,j] = signal / noise
      end
    end
  end
  SNR[:,:,:] .= mean(SNR,dims=3)
  return
end

function calculateSNRCustomSF(f::BrukerFile,fgFrames::Array,bgFramesFull::Array,bgFramesHalf::Array)
  SNR = zeros(rxNumFrequencies(f),rxNumChannels(f),1)
  for j=1:1
    for r=1:rxNumChannels(f)
      for k=1:rxNumFrequencies(f)
        meanBG = mean(bgFramesFull[:,k,r,j])
        signal = maximum(abs.(fgFrames[:,k,r,j].-meanBG))
	meanBGHalf = mean(bgFramesHalf[:,k,r,j])
        noise = sqrt(var(bgFramesFull[:,k,r,j].-meanBGHalf))#[:,k,r,j]))
        SNR[k,r,j] = signal / noise
      end
    end
  end
  SNR[:,:,:] .= mean(SNR,dims=3)
  return SNR
end

function converttoreal(S::AbstractArray{Complex{T},2}) where {T}
  N = size(S,1)
  M = size(S,2)
  S = reshape(reinterpret(T,vec(S)),(2*N,M))
  for l=1:M
    tmp = S[:,l]
    S[1:N,l] = tmp[1:2:end]
    S[N+1:end,l] = tmp[2:2:end]
  end
  return reshape(S,(N,2*M))
end
