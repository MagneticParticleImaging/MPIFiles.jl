function getInterpolator(A::AbstractArray{T,3}, grid::RegularGridPositions) where T

  tmp = ( size(A,1)==1 ? NoInterp() : BSpline(Linear()),
          size(A,2)==1 ? NoInterp() : BSpline(Linear()),
          size(A,3)==1 ? NoInterp() : BSpline(Linear()) )

  itp = extrapolate(interpolate(A, tmp), 0.0)
  sitp = scale(itp, size(A,1)==1 ? Base.OneTo(1) : range(grid,1), size(A,2)==1 ? Base.OneTo(1) : range(grid,2), size(A,3)==1 ? Base.OneTo(1) : range(grid,3))
  return sitp
end

function interpolate(A::AbstractArray{T,3}, origin::RegularGridPositions,
                        target::RegularGridPositions) where T

  sitp = getInterpolator(A, origin)
  N = target.shape
  AInterp = zeros(eltype(A), N[1], N[2], N[3])
  # We need to handle origin grids with size 1 in one dimension differently, as the interpolation object can not be scaled correctly to use positions instead of indices.
  # Therefore the only valid interpolation position is with index 1. We search for positions close enough to the origin and use 1 there
  ranges = [if origin.shape[i]==1
        tmp = zeros(N[i])
        tmp[isapprox.(range(target,i),origin.center[i], atol=1e-12)] .= 1
        tmp
      else
         range(target,i)
      end for i in 1:3]
  return _interpolate_inner(AInterp,N,sitp,ranges...)
end

function _interpolate_inner(AInterp,N,sitp,rx,ry,rz)
  for nz=1:N[3]
    for ny=1:N[2]
      for nx=1:N[1]
          AInterp[nx,ny,nz] = sitp(rx[nx],ry[ny],rz[nz])
      end
    end
  end

  return AInterp
end
