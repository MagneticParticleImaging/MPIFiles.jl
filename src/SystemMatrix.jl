export getSystemMatrix, getSystemMatrixReshaped, calculateSystemMatrixSNR

"""
  getSystemMatrix(f, [neglectBGFrames]; kargs...) => Array{ComplexF32,4}

Load the system matrix in frequency domain

Supported keyword arguments:
* frequencies
* bgCorrection
* loadasreal
"""
function getSystemMatrix(f::MPIFile, frequencies=nothing;
                         bgCorrection=false, loadasreal=false,
                         tfCorrection=rxHasTransferFunction(f), 
                         numPeriodAverages=1, numPeriodGrouping=1, kargs...)

  if isnothing(frequencies)
    frequencies = filterFrequencies(f)
  end

  if measIsFastFrameAxis(f) && measIsFourierTransformed(f)
    # This is the fast path if the SM lays on disk in the required format
    data = systemMatrix(f, frequencies, bgCorrection)
  else
    data = reshape(getMeasurementsFD(f, frequencies=frequencies, sortFrames=true, bgCorrection=bgCorrection,
              spectralLeakageCorrection=false, transposed=true, tfCorrection=false,
              numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping), Val(2))
  end
  
  S = map(ComplexF32, data)

  if tfCorrection && !measIsTFCorrected(f)
    tf = rxTransferFunction(f)
    if !isnothing(tf)
      _corrTFSF(S,tf[rowsToSubsampledRows(f,frequencies)])
      @info "System matrix has been corrected with a Transfer Function. Name of TF: $(rxTransferFunctionFileName(f))"
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

function calculateSystemMatrixSNR(f::MPIFile, S::Array; numPeriodAverages=1, numPeriodGrouping=1)
  J = div(acqNumPeriodsPerFrame(f),numPeriodAverages*numPeriodGrouping)
  R = rxNumChannels(f)
  K = rxNumFrequencies(f,numPeriodGrouping)
  N = acqNumFGFrames(f)

  SNR = zeros(K, R, J)

  gridSize = tuple(calibSize(f)...)

  calculateSystemMatrixSNRInner(S, SNR, J, R, K, N, gridSize)
  return SNR
end

function calculateSystemMatrixSNRInner(S, SNR, J, R, K, N, gridSize::NTuple{D,Int}) where D
  for j=1:J
    for r=1:R
      # precalc on entire FOV
      for k=1:K
        SBG = S[(N+1):end,k,r,j]
        SFG = S[1:N,k,r,j]
        meanBG = mean(SBG)
        noise = mean(abs.(SBG .- meanBG))
        signal = median( abs.(SFG.-meanBG) ) 
        SNR[k,r,j] = signal / noise
      end
      # generate mask representing signal region 
      idx = sortperm(vec(SNR[:,r,j]), rev=true)
      mask = zeros(Bool, N)
      for q = 1:20
        SFG = abs.(S[1:N,idx[q],r,j])
        maxSFG = maximum(SFG)
        mask[ SFG .> maxSFG*0.90 ] .= true
      end
      #@info sum(mask)/length(mask)
      # calc SNR on mask
      for k=1:K
        SBG = S[(N+1):end,k,r,j]
        SFG = S[1:N,k,r,j]
        meanBG = mean(SBG)
        noise = mean(abs.(SBG .- meanBG))
        κ = abs.(SFG.-meanBG) 
        #κ_ = mapwindow(median!, reshape(κ, gridSize), ntuple(d->5,D))
        #signal = maximum( κ_ ) 
        signal = median( κ[mask .== true] ) 
        #signal = maximum( κ[mask .== true] ) 
        if signal > 3.5*noise
          phaseMaskA = angle.(SFG.-meanBG) .> 0
          phaseMaskB = angle.(SFG.-meanBG) .<= 0
	  phaseMask = sum(phaseMaskA) > sum(phaseMaskB) ? phaseMaskA : phaseMaskB
	  #signal = maximum( κ[mask .== true] )
          κmax = maximum(κ)
          signal = median( κ[(κ .> κmax*0.9 ) ] ) 
          #signal = mean( κ[ mask .&&  phaseMask ] ) 
        end
        SNR[k,r,j] = signal / noise
      end
    end
  end
  SNR[:,:,:] .= median(SNR,dims=3)
  return
end
#=
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
=#
function calculateSNRCustomSF(f::BrukerFile,fgFrames::Array,bgFramesFull::Array,bgFramesHalf::Array)
  SNR = zeros(rxNumFrequencies(f),rxNumChannels(f),1)
  for j=1:1
    for r=1:rxNumChannels(f)
      for k=1:rxNumFrequencies(f)
        diffBG = diff(bgFramesFull[:,k,r,j])
        meanBG = mean(bgFramesFull[:,k,r,j])
        signal = maximum(abs.(fgFrames[:,k,r,j].-meanBG))
        noise = mean(abs.(diffBG))
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
