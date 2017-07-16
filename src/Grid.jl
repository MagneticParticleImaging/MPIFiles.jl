using Unitful

import Base: getindex, length, convert

export AbstractGrid, RegularGrid, MeanderingGrid, ArbitraryGrid

abstract AbstractGrid

type RegularGrid{T} <: AbstractGrid
  shape::Vector{Int}
  fov::Vector{T}
  center::Vector{T}
  #RegularGrid(shape,fov,center=[0.0,0.0,0.0]) = new(shape,fov,center)
end

getindex(grid::RegularGrid, i::Integer) =
  (collect(ind2sub(tuple(grid.shape...), i)).-0.5.-grid.shape./2)./
              grid.shape.*grid.fov .- grid.center


type MeanderingGrid{T} <: AbstractGrid
  shape::Vector{Int}
  fov::Vector{T}
  center::Vector{T}
end

function getindex(grid::MeanderingGrid, i::Integer)
  idx = collect(ind2sub(tuple(grid.shape...), i))
  for d=2:3
    if iseven(idx[d])
      idx[d-1] = grid.shape[d-1] - idx[d-1] + 1
    end
  end

  return (idx.-0.5.-grid.shape./2)./
              grid.shape.*grid.fov .- grid.center
end

length(grid::Union{RegularGrid,MeanderingGrid}) = prod(grid.shape)

type ArbitraryGrid{T} <: AbstractGrid
  positions::Matrix{T}
  #RegularGrid(shape,fov,center=[0.0,0.0,0.0]) = new(shape,fov,center)
end

getindex(grid::ArbitraryGrid, i::Integer) = grid.positions[:,i]

function convert(::Type{ArbitraryGrid}, grid::AbstractGrid)
  positions = zeros(3,length(grid))
  for i=1:length(grid)
    positions[:,i] = grid[i]
  end
  return ArbitraryGrid(positions)
end








### Old interface ###


export createRandPositions, calcNumberOfFrames
export createCartesianSF

const timeSpanPerFrame3D = (21.542*10.0^-3)u"s" #seconds
const overlapFrames = 30

const backGroundPos = [220 0 0];

@doc "Creates n random positions in x,y,z, in mm waiting time will be calculated from frames,
optional seed for MersenneTwister generator"->
function createRandPositions(FOV::Array{typeof(1.0u"mm"),1}, numPos::Unsigned, numFramesToAcquire::Unsigned;
    offset::Array{typeof(1.0u"mm"),1} = [0.0u"mm"; 0.0u"mm"; 0.0u"mm"], waitLingerTime::typeof(1.0u"s") = 0.0u"s", seed::Unsigned = UInt64(0))
  if length(FOV) != 3 || length(offset) != 3
    error("FOV must have size 3 for random positions in x,y,z")
  end
  mersenneTwister = MersenneTwister(rand(UInt64,1)[1])
  if seed != UInt64(0)
    mersenneTwister = MersenneTwister(seed)
  end
  # scale to FOV
  randPos = Array{typeof(1.0u"mm"),2}(numPos, 3);
  rP = rand(mersenneTwister, numPos, 3)
  randPos[:,1]=rP[:,1].*FOV[1]
  randPos[:,2]=rP[:,2].*FOV[2]
  randPos[:,3]=rP[:,3].*FOV[3]

  # move to center of FOV
  randPos[:,1]=randPos[:,1]-(FOV[1]/2) + offset[1]
  randPos[:,2]=randPos[:,2]-(FOV[2]/2) + offset[2]
  randPos[:,3]=randPos[:,3]-(FOV[3]/2) + offset[3]

  # set time to wait 2x15=30 frames overlap
  randPosTime = ones(numPos) .* numFramesToAcquire * timeSpanPerFrame3D +
                  overlapFrames * timeSpanPerFrame3D + waitLingerTime;

  framesTotal = calcNumberOfFrames(numPos, waitLingerTime, framesPerPosition = numFramesToAcquire)
  println("Est. required frames in BrukerSoftware: $(framesTotal)")
  randPos, randPosTime, numFramesToAcquire, seed
end

function calcNumberOfFrames(numPos, waitingTime::typeof(1.0u"s"); framesPerPosition=100)
  latency = 3000
  movementTimePerPosition = 2000
  numframes = (latency + movementTimePerPosition* numPos + Float64(ustrip(uconvert(u"ms",timeSpanPerFrame3D)))
  * framesPerPosition * numPos + Float64(ustrip(uconvert(u"ms", waitingTime))) * numPos) / Float64(ustrip(uconvert(u"ms",timeSpanPerFrame3D)))
  if numframes > 100000
    warn("Number of frames exceed maximal setting of frames 100000")
  end
  numframes
end

"""Creates cartesian grid in x,y,z, in mm with `gridsize`, `backGroundInc` and `measureTime`"""
function createCartesianSF(gridSize, backGroundInc, measureTime)
  numPoints = prod(gridSize) + prod(gridSize[2:3]) + 1
  coords = Array{Float64,2}(numPoints,3)
  i=1
  coords[i,:] = backGroundPos
  i+=1
  for k=1:gridSize[3]
    for l=1:gridSize[2]
      for m=1:gridSize[1]
        coords[i,:] = [m l k].-[gridSize[1]/2 gridSize[2]/2 gridSize[3]/2]
        i+=1
        if mod(m,backGroundInc)==zero(0)
          coords[i,:] = backGroundPos
          i+=1
        end
      end
    end
  end
  timeSpan = ones(numPoints) * measureTime
  coords = coords *u"mm"
  return coords, timeSpan
end
