export getSystemMatrix

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

function getSystemMatrix(f::MPIFile, frequencies; bgcorrection=false, loadasreal=false,
                         loadas32bit=false, kargs...)
  data = procData(f, frequencies)

  S = loadas32bit ? map(Complex64, data) : map(Complex128, data)

  if loadasreal
    return converttoreal(S)
  else
    return S
  end
end

# TODO: Sparse system matrix loading
