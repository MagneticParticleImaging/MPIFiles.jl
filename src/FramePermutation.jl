# The following function allow to permute the frames
export meanderingFramePermutation, fullFramePermutation

# TODO the following requires a test
function meanderingFramePermutation(f::MPIFile)
  N = tuple(calibSize(f)...)

  perm = Array{Int}(N)
  for i in CartesianRange(N)
    idx = [i[k] for k=1:length(i)]
    for d=2:3
      if isodd(sum(idx[d:3])-length(idx[d:3]))
        idx[d-1] = N[d-1] + 1 - idx[d-1]
      end
    end
    perm[i] = sub2ind(N,idx...)
  end
  return vec(perm)
end

function fullFramePermutation(f::MPIFile, meandering::Bool)
  perm1=cat(1,measFGFrameIdx(f),measBGFrameIdx(f))
  if meandering
    perm2=cat(1,meanderingFramePermutation(f),
               (length(perm1)-acqNumBGFrames(f)+1):length(perm1))
    permJoint = perm1[perm2]
    return permJoint
  else
    return perm1
  end
end
