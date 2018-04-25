export getSystemMatrix, getSystemMatrixReshaped

function converttoreal{T}(S::AbstractArray{Complex{T},2})
  N = size(S,1)
  M = size(S,2)
  S = reinterpret(T,S,(2*N,M))
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
  S = map(Complex64, data)

  if loadasreal
    return converttoreal(S)
  else
    return S
  end
end

function getSystemMatrixReshaped(f::MPIFile; kargs...)
  return reshape(getSystemMatrix(f),:,rxNumFrequencies(f),rxNumChannels(f),acqNumPeriodsPerFrame(f))
end
