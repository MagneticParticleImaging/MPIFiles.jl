export Positions, GridPositions, NestedPositions, RegularGridPositions, ChebyshevGridPositions,
       MeanderingGridPositions, UniformRandomPositions, ArbitraryPositions, SortedPositions,
       SphericalTDesign, BreakpointPositions, BreakpointGridPositions
export SpatialDomain, AxisAlignedBox, Ball
export loadTDesign, getPermutation
export fieldOfView, fieldOfViewCenter, shape
export posToIdx, posToLinIdx, spacing, isSubgrid, deriveSubgrid, toDict
export axesToRegularGridPositions

abstract type Positions{T, D} end
abstract type GridPositions{T, D} <: Positions{T, D} end
ndims(pos::Positions{T, D}) where {T, D} = D

abstract type NestedPositions{T, D, G <: Positions{T, D}} <: Positions{T, D} end
parent(position::NestedPositions) = throw(error("$(typeof(position)) must implement `parent`"))

function Positions(file::HDF5.File)
  if haskey(file, "/positionsBreakpoint")
    return BreakpointPositions(file)
  end

  typ = read(file, "/positionsType")
  if typ == "RegularGridPositions"
    positions = RegularGridPositions(file)
  elseif typ == "ChebyshevGridPositions"
    positions = ChebyshevGridPositions(file)
  elseif typ == "SphericalTDesign"
    positions = SphericalTDesign(file)
  elseif typ == "UniformRandomPositions"
    positions = UniformRandomPositions(file)
  elseif typ == "ArbitraryPositions"
    positions = ArbitraryPositions(file)
  else
    throw(ErrorException("No grid found to load from $file"))
  end

  if haskey(file, "/positionsMeandering") && typ in ["RegularGridPositions","ChebyshevGridPositions"] && read(file, "/positionsMeandering") == Int8(1)
    positions = MeanderingGridPositions(positions)
  end

  return positions
end

function Positions(params::Dict)
  type = params["type"]
  if type == "RegularGridPositions"
    positions = RegularGridPositions(params)
  elseif type == "MeanderingGridPositions"
    positions = MeanderingGridPositions(params)
  elseif type == "ChebyshevGridPositions"
    positions = ChebyshevGridPositions(params)
  elseif type == "BreakpointPositions"
    positions = BreakpointPositions(params)
  elseif type == "SphericalTDesign"
    positions = SphericalTDesign(params)
  elseif type == "UniformRandomPositions"
    positions = UniformRandomPositions(params)
  elseif type == "ArbitraryPositions"
    positions = ArbitraryPositions(params)
  elseif type == "SortedPositions"
    positions = SortedPositions(params)
  elseif type == "TubularRegularGridPositions"
    positions = TubularRegularGridPositions(params)
  elseif type == "SubsampledPositions"
    positions = SubsampledPositions(params)
  else
    throw(ErrorException("No grid found to load from dict $params"))
  end

  return positions
end

# Cartesian grid
struct RegularGridPositions{T, D} <: GridPositions{T, D}
  shape::SVector{D, Int64}
  fov::SVector{D, T}
  center::SVector{D, T}
  sign::SVector{D, Int64}
end
RegularGridPositions(shape, fov, center, signs) = RegularGridPositions(SVector{length(shape)}(shape), SVector{length(fov)}(fov), SVector{length(center)}(center), SVector{length(signs)}(signs))

function range(grid::RegularGridPositions, dim::Int)
  if grid.shape[dim] > 1
    sp = spacing(grid)
    return range(grid.center[dim] - grid.sign[dim]*(grid.fov[dim]/2 - sp[dim]/2),
                 step=grid.sign[dim]*sp[dim], length=grid.shape[dim])
  else
    return range(start=grid.center[dim],stop=grid.center[dim],step=oneunit(grid.center[dim]))
  end
end
Base.axes(grid::RegularGridPositions) = tuple([range(grid, i) for i in 1:ndims(grid)]...)

RegularGridPositions(shape, fov, center) = RegularGridPositions(shape, fov, center, ones(Int,length(shape)))

function RegularGridPositions(file::HDF5.File)
  shape = read(file, "/positionsShape")
  fov = read(file, "/positionsFov")*Unitful.m
  center = read(file, "/positionsCenter")*Unitful.m

  return RegularGridPositions(shape, fov, center)
end

function RegularGridPositions(params::Dict)
  if params["type"] != "RegularGridPositions"
    throw(ArgumentError("Unexpected positions type $(params["type"])"))
  end

  shape = params["shape"]

  fov = params["fov"]
  if !(eltype(fov) <: Quantity)
    fov = fov.*Unitful.m
  end

  center = params["center"]
  if !(eltype(center) <: Quantity)
    center = center.*Unitful.m
  end

  return RegularGridPositions(shape, fov, center)
end

# Find a joint grid
function RegularGridPositions(positions::Vector{T}) where T<:RegularGridPositions
  posMin = positions[1].center .- 0.5*positions[1].fov
  posMax = positions[1].center .+ 0.5*positions[1].fov
  minSpacing = spacing(positions[1])
  for position in positions
    sp = spacing(position)
    for d=1:length(posMin)
      posMin[d] = min(posMin[d], position.center[d] - 0.5*position.fov[d])
      posMax[d] = max(posMax[d], position.center[d] + 0.5*position.fov[d])
      minSpacing[d] = min(minSpacing[d],sp[d])
    end
  end
  center = (posMin .+ posMax)/2
  fov = posMax .- posMin
  shape = round.(Int64,fov./minSpacing)
  fov = shape .* minSpacing
  return RegularGridPositions(shape, fov, center)
end

function isSubgrid(grid::RegularGridPositions, subgrid::RegularGridPositions)
  if any(fieldOfView(grid) .- fieldOfView(subgrid) .< 0) ||
     any(spacing(grid) .!= spacing(subgrid))
    return false
  else
    centerPosIdx = posToIdxFloat(grid,subgrid[ones(Int,length(subgrid.shape))])
    return all(isapprox.(centerPosIdx, round.(Int,centerPosIdx) ,rtol=1e-5))
  end
end

function deriveSubgrid(grid::RegularGridPositions, subgrid::RegularGridPositions)
  minI = ones(Int,length(subgrid.shape))
  maxI = copy(subgrid.shape)
  for d=1:length(minI)
    if subgrid.sign[d] == -1
      minI[d] = subgrid.shape[d]-minI[d]+1
      maxI[d] = subgrid.shape[d]-maxI[d]+1
    end
  end
  minPos = subgrid[ minI ]
  maxPos = subgrid[ maxI ]

  minIdx = posToIdx(grid,minPos)
  maxIdx = posToIdx(grid,maxPos)
  #shp = maxIdx-minIdx+ones(Int,length(subgrid.shape))
  shp = shape(subgrid)
  #center = (grid[minIdx].+grid[maxIdx])/2
  # TODO round properly
  center = (grid[minIdx].+grid[minIdx.+shp.-1])/2
  fov = shp.*spacing(grid)
  return RegularGridPositions(shp,fov,center,subgrid.sign)
end

function write(file::HDF5.File, positions::RegularGridPositions)
  write(file,"/positionsType", "RegularGridPositions")
  write(file, "/positionsShape", Array(positions.shape))
  write(file, "/positionsFov", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.fov)))) )
  write(file, "/positionsCenter", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center)))) )
end

function toDict(positions::RegularGridPositions)
  params = Dict{String,Any}()
  params["type"] = "RegularGridPositions"
  params["shape"] = positions.shape
  params["fov"] = Float64.(ustrip.(uconvert.(Unitful.m, positions.fov)))
  params["center"] = Float64.(ustrip.(uconvert.(Unitful.m, positions.center)))
  return params
end

function getindex(grid::RegularGridPositions, i::Integer)
  if i>length(grid) || i<1
     throw(BoundsError(grid,i))
  end

  #idx = collect(ind2sub(tuple(shape(grid)...), i))
  if length(grid.shape) == 1 #Very ugly but improves compile time
    idx = [i]
  elseif length(grid.shape) == 2
    idx = collect(Tuple((CartesianIndices(tuple(grid.shape[1],grid.shape[2])))[i]))
  elseif length(grid.shape) == 3
    idx = collect(Tuple((CartesianIndices(tuple(grid.shape[1],grid.shape[2],grid.shape[3])))[i]))
  else
    idx = collect(Tuple((CartesianIndices(tuple(grid.shape...)))[i]))
  end

  for d=1:length(idx)
    if grid.sign[d] == -1
      idx[d] = grid.shape[d]-idx[d]+1
    end
  end
  return ((-shape(grid).+(2 .*idx.-1))./shape(grid)).*fieldOfView(grid)./2 + fieldOfViewCenter(grid)
end

function getindex(grid::RegularGridPositions, idx::Vector{T}) where T<:Number
  for d=1:length(idx)
    if grid.sign[d] == -1
      idx[d] = grid.shape[d]-idx[d]+1
    end
  end
  return 0.5.*fieldOfView(grid) .* (-1 .+ (2 .* idx .- 1) ./ shape(grid)) .+ fieldOfViewCenter(grid)
end

function getindex(grid::RegularGridPositions, idx::CartesianIndex)
  return getindex(grid, LinearIndices(tuple(grid.shape...))[idx])
end

function posToIdxFloat(grid::RegularGridPositions, pos)
  idx = 0.5 .* (shape(grid) .* ((pos .- fieldOfViewCenter(grid)) ./
              ( 0.5 .* fieldOfView(grid) ) .+ 1) .+ 1)
  idx = [isnan(val) ? one(eltype(idx)) : val for val in idx]
  return idx
end

function posToIdx(grid::RegularGridPositions, pos)
  idx = round.(Int64, posToIdxFloat(grid,pos))
  for d=1:length(idx)
    if grid.sign[d] == -1
      idx[d] = grid.shape[d]-idx[d]+1
    end
  end
  return idx
end

function posToLinIdx(grid::RegularGridPositions, pos)
  return (LinearIndices(tuple(shape(grid)...)))[posToIdx(grid,pos)...]
end

"""
Create a `RegularGridPositions` object from three axes definition the coordinates of the positions.

Automatically checks that the axes all have regular spacing.

- `axesToRegularGridPositions(axs::Tuple)`
- `axesToRegularGridPositions(x, y, z)`
"""
axesToRegularGridPositions(axs::Tuple) = axesToRegularGridPositions(axs...)
function axesToRegularGridPositions(x, y, z)
  allsame(ax) = length(ax)==0 || all(≈(first(ax)),ax)
  if !allsame(diff(x)) || !allsame(diff(y)) || !allsame(diff(z)) 
    error("Axes do not produce a regular grid!")
  end
  shape = [length(ax) for ax in [x,y,z]]
  fov = [abs(ax[end]-ax[1])*length(ax)/max(1,length(ax)-1) for ax in [x,y,z]]
  center = [mean(ax) for ax in [x,y,z]]
  sign = [if ax[end]<ax[1]; -1 else 1 end for ax in [x,y,z]]
  return RegularGridPositions(shape, fov, center, sign)
end


# Chebyshev Grid
struct ChebyshevGridPositions{T, D, S} <: GridPositions{T, D}
  shape::SVector{D, Int64}
  fov::SVector{D, S}
  center::SVector{D, T}
end
ChebyshevGridPositions(shape, fov, center) = ChebyshevGridPositions(SVector{length(shape)}(shape), SVector{length(fov)}(fov), SVector{length(center)}(center))

function write(file::HDF5.File, positions::ChebyshevGridPositions)
  write(file,"/positionsType", "ChebyshevGridPositions")
  write(file, "/positionsShape", Array(positions.shape))
  write(file, "/positionsFov", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.fov)))) )
  write(file, "/positionsCenter", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center)))) )
end

function toDict(positions::ChebyshevGridPositions)
  params = Dict{String,Any}()
  params["type"] = "ChebyshevGridPositions"
  params["shape"] = Array(positions.shape)
  params["fov"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.fov))))
  params["center"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center))))
  return params
end

function ChebyshevGridPositions(file::HDF5.File)
  shape = read(file, "/positionsShape")
  fov = read(file, "/positionsFov")*Unitful.m
  center = read(file, "/positionsCenter")*Unitful.m
  return ChebyshevGridPositions(shape,fov,center)
end

function ChebyshevGridPositions(params::Dict)
  if params["type"] != "ChebyshevGridPositions"
    throw(ArgumentError("Unexpected positions type $(params["type"])"))
  end

  shape = params["shape"]

  fov = params["fov"]
  if !(eltype(fov) <: Quantity)
    fov = fov.*Unitful.m
  end

  center = params["center"]
  if !(eltype(center) <: Quantity)
    center = center.*Unitful.m
  end

  return ChebyshevGridPositions(shape,fov,center)
end

function getindex(grid::ChebyshevGridPositions, i::Integer)
  if i>length(grid) || i<1
    throw(BoundsError(grid,i))
  else
    idx = collect(Tuple(CartesianIndices(tuple(shape(grid)...))[i]))
    return -cos.((idx .- 0.5) .* pi ./ shape(grid)) .* fieldOfView(grid) ./ 2 .+ fieldOfViewCenter(grid)
  end
end

# Meander regular grid positions
struct MeanderingGridPositions{T, D, G <: GridPositions{T, D}} <: NestedPositions{T, D, G}
  grid::G
end
parent(grid::MeanderingGridPositions) = grid.grid

function MeanderingGridPositions(file::HDF5.File)
  typ = read(file, "/positionsType")
  if typ == "RegularGridPositions"
    grid = RegularGridPositions(file)
    return MeanderingGridPositions(grid)
  elseif typ == "ChebyshevGridPositions"
    grid = ChebyshevGridPositions(file)
    return MeanderingGridPositions(grid)
  end
end

function MeanderingGridPositions(params::Dict)
  pos = Positions(params["positions"])
  return MeanderingGridPositions(pos)
end

function write(file::HDF5.File, positions::MeanderingGridPositions)
  write(file,"/positionsMeandering", Int8(1))
  write(file, positions.grid)
end

function toDict(positions::MeanderingGridPositions)
  params = Dict{String,Any}()
  params["positions"] = toDict(positions.grid)
  params["type"] = "MeanderingGridPositions"
  return params
end

function indexPermutation(grid::MeanderingGridPositions, i::Integer)
  dims = tuple(shape(grid)...)
  ndims = length(dims)
  idx = collect(Tuple(CartesianIndices(dims)[i]))
    for d=1:ndims-1
      if isodd(sum(idx[d+1:end])-(ndims-d))
      idx[d] = -idx[d] + shape(grid)[d] + 1
    end
  end
  linidx = (LinearIndices(dims))[idx...]
end

function getindex(grid::MeanderingGridPositions, i::Integer)
  iperm = indexPermutation(grid,i)
  return grid.grid[iperm]
end

function getPermutation(grid::MeanderingGridPositions)
  N = length(grid)
  perm = Array{Int}(undef,N)

  for i in eachindex(perm)
    perm[i] = indexPermutation(grid,i)
  end
  return vec(perm)
end

struct BreakpointPositions{T, D, G} <: NestedPositions{T, D, G}
  grid::G
  breakpointIndices::Vector{Int64}
  breakpointPosition::SVector{D, T}
end
const BreakpointGridPositions = BreakpointPositions
BreakpointPositions(grid, indices, pos) = BreakpointPositions(grid, indices, SVector{length(pos)}(pos))
parent(grid::BreakpointPositions) = grid.grid

function BreakpointPositions(file::HDF5.File)
  typ = read(file, "/positionsType")
  breakpointPosition = read(file, "/positionsBreakpoint") * Unitful.m
  breakpointIndices = read(file, "/indicesBreakpoint")

  if typ == "MeanderingGridPositions"
    grid = MeanderingGridPositions(file)
    return BreakpointPositions(grid,breakpointIndices, breakpointPosition)
  elseif typ == "RegularGridPositions"
    grid = RegularGridPositions(file)
    return BreakpointPositions(grid, breakpointIndices, breakpointPosition)
  elseif typ == "ChebyshevGridPositions"
    grid = ChebyshevGridPositions(file)
    return BreakpointPositions(grid, breakpointIndices, breakpointPosition)
  end
end

function BreakpointPositions(params::Dict)
  breakpointPosition = params["breakpoint"]
  if !(eltype(breakpointPosition) <: Quantity)
    breakpointPosition = breakpointPosition.*Unitful.m
  end

  breakpointIndices = params["indices"]
  positions = Positions(params["positions"])
  return BreakpointGridPositions(positions, breakpointIndices, breakpointPosition)
end

function write(file::HDF5.File, positions::BreakpointPositions)
  write(file,"/positionsBreakpoint", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.breakpointPosition)))))
  write(file,"/indicesBreakpoint", positions.breakpointIndices)
  write(file, positions.grid)
end

function toDict(positions::BreakpointPositions)
  params = Dict{String, Any}()
  params["positions"] = toDict(positions.grid)
  params["type"] = "BreakpointPositions"
  params["breakpoint"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.breakpointPosition))))
  params["indices"] = positions.breakpointIndices
  return params
end

function getmask(grid::BreakpointPositions)
  bgind=grid.breakpointIndices
  mask = zeros(Bool, length(grid.grid)+length(bgind))
  mask[bgind] .= true
  return mask
end

function getindex(grid::BreakpointPositions, i::Integer)

  bgind=grid.breakpointIndices

  if i>(length(grid.grid)+length(bgind)) || i<1
    return throw(BoundsError(grid,i))
  elseif any(i .== bgind)
    return grid.breakpointPosition
  else
    pastBgind = sum(i .> bgind)
    return grid.grid[i-pastBgind]
  end
end

# Uniform random distributed positions
abstract type SpatialDomain{S, T} end

struct AxisAlignedBox{S, T} <: SpatialDomain{S, T}
  fov::Vector{S}
  center::Vector{T}
end

function write(file::HDF5.File, domain::AxisAlignedBox)
  write(file, "/positionsDomain", "AxisAlignedBox")
  write(file, "/positionsDomainFieldOfView", Float64.(ustrip.(uconvert.(Unitful.m, domain.fov))) )
  write(file, "/positionsDomainCenter", Float64.(ustrip.(uconvert.(Unitful.m, domain.center))) )
end

function AxisAlignedBox(file::HDF5.File)
  fov = read(file, "/positionsDomainFieldOfView")*Unitful.m
  center = read(file, "/positionsDomainCenter")*Unitful.m
  return AxisAlignedBox(fov,center)
end

struct Ball{S, T} <: SpatialDomain{S, T}
  radius::S
  center::Vector{T}
end

function write(file::HDF5.File, domain::Ball)
  write(file, "/positionsDomain", "Ball")
  write(file, "/positionsDomainRadius", Float64.(ustrip.(uconvert.(Unitful.m, domain.radius))) )
  write(file, "/positionsDomainCenter", Float64.(ustrip.(uconvert.(Unitful.m, domain.center))) )
end

function Ball(file::HDF5.File)
  radius = read(file, "/positionsDomainRadius")*Unitful.m
  center = read(file, "/positionsDomainCenter")*Unitful.m
  return Ball(radius,center)
end


mutable struct UniformRandomPositions{T, D, S <: SpatialDomain} <: Positions{T, D}
  N::UInt
  seed::UInt32
  domain::S
end
UniformRandomPositions(N, seed, domain::SpatialDomain{S, T}) where {S,T} = UniformRandomPositions{T, 3, typeof(domain)}(N, seed, domain)

radius(rpos::UniformRandomPositions{T, D, <:Ball}) where {T, D} = rpos.domain.radius
seed(rpos::UniformRandomPositions) = rpos.seed

function getindex(rpos::UniformRandomPositions{T, D, <:AxisAlignedBox}, i::Integer) where {T, D}
  if i>length(rpos) || i<1
    throw(BoundsError(rpos,i))
  else
    # make sure Positions are randomly generated from given seed
    rng = StableRNG(seed(rpos))
    rP = rand(rng, 3, i)[:,i]
    return (rP.-0.5).*fieldOfView(rpos)+fieldOfViewCenter(rpos)
  end
end

function getindex(rpos::UniformRandomPositions{T, Dim, <:Ball}, i::Integer) where {T, Dim}
  if i>length(rpos) || i<1
    throw(BoundsError(rpos,i))
  else
    # make sure Positions are randomly generated from given seed
    rng = StableRNG(seed(rpos))
    D = rand(rng, i)[i]
    P = randn(rng, 3, i)[:,i]
    return radius(rpos)*D^(1/3)*normalize(P)+fieldOfViewCenter(rpos)
  end
end

function write(file::HDF5.File, positions::UniformRandomPositions)
  write(file, "/positionsType", "UniformRandomPositions")
  write(file, "/positionsN", positions.N)
  write(file, "/positionsSeed", positions.seed)
  write(file, positions.domain)
end

function UniformRandomPositions(file::HDF5.File)
  N = read(file, "/positionsN")
  seed = read(file, "/positionsSeed")
  dom = read(file,"/positionsDomain")
  if dom=="Ball"
    domain = Ball(file)
    return UniformRandomPositions(N,seed,domain)
  elseif dom=="AxisAlignedBox"
    domain = AxisAlignedBox(file)
    return UniformRandomPositions(N,seed,domain)
  else
    throw(ErrorException("No method to read domain $domain"))
  end
end

# TODO fix conversion methods
#=
function convert(::Type{UniformRandomPositions}, N::Integer,seed::UInt32,fov::Vector{S},center::Vector{T}) where {S,T<:Unitful.Length}
  if N<1
    throw(DomainError())
  else
    uN = convert(UInt,N)
    return UniformRandomPositions(uN,seed,fov,center)
  end
end

function convert(::Type{UniformRandomPositions}, N::Integer,fov::Vector,center::Vector)
  return UniformRandomPositions(N,rand(UInt32),fov,center)
end
=#

# Tubular cartesian grid
export TubularRegularGridPositions
struct TubularRegularGridPositions{T, D} <: GridPositions{T, D}
  shape::SVector{D, Int64}
  fov::SVector{D, T}
  center::SVector{D, T}
  "Main axis of the tube; only effective in 3D grids"
  mainAxis::Int64
  "Radius-defining axis of the tube"
  radiusAxis::Int64
end
function TubularRegularGridPositions(shape, fov, center, mainAxis, radius) 
  TubularRegularGridPositions(SVector{length(shape)}(shape), SVector{length(fov)}(fov), SVector{length(center)}(center), mainAxis, radius)
end

function TubularRegularGridPositions(params::Dict)
  shape = params["shape"]

  fov = params["fov"]
  if !(eltype(fov) <: Quantity)
    fov = fov.*Unitful.m
  end

  center = params["center"]
  if !(eltype(center) <: Quantity)
    center = center.*Unitful.m
  end

  mainAxis = params["mainAxis"]
  radiusAxis = params["radiusAxis"]

  return TubularRegularGridPositions(shape, fov, center, mainAxis, radiusAxis)
end

function TubularRegularGridPositions(file::HDF5.File)
  shape = read(file, "/positionsShape")
  fov = read(file, "/positionsFov")*Unitful.m
  center = read(file, "/positionsCenter")*Unitful.m
  mainAxis = read(file, "positionsMainAxis")
  radiusAxis = read(file, "positionsRadiusAxis")

  return TubularRegularGridPositions(shape, fov, center, mainAxis, radiusAxis)
end

length(grid::TubularRegularGridPositions) = length(filteredPositions(grid))

export radius
radius(grid::TubularRegularGridPositions) = grid.fov[grid.radiusAxis] / 2

function write(file::HDF5.File, positions::TubularRegularGridPositions)
  write(file, "/positionsType", "TubularRegularGridPositions")
  write(file, "/positionsShape", Array(positions.shape))
  write(file, "/positionsFov", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.fov)))) )
  write(file, "/positionsCenter", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center)))) )
  write(file, "/positionsMainAxis", positions.mainAxis)
  write(file, "/positionsRadiusAxis", positions.radiusAxis)
end

function toDict(positions::TubularRegularGridPositions)
  params = Dict{String,Any}()
  params["type"] = "TubularRegularGridPositions"
  params["shape"] = Array(positions.shape)
  params["fov"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.fov))))
  params["center"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center))))
  params["mainAxis"] = positions.mainAxis
  params["radiusAxis"] = positions.radiusAxis
  return params
end

function filteredPositions(grid::TubularRegularGridPositions)
  if length(grid.shape) == 1 #Very ugly but improves compile time
    cartIndices = CartesianIndices(tuple(grid.shape[1]))
    return [Tuple(idx) for idx in cartIndices] # Applying a radius in 1D is not possible
  elseif length(grid.shape) == 2
    cartIndices = CartesianIndices(tuple(grid.shape[1], grid.shape[2]))
    return [Tuple(idx) for idx in cartIndices if norm(collect(Tuple(idx)) .- grid.shape./2)]
  elseif length(grid.shape) == 3
    cartIndices = CartesianIndices(tuple(grid.shape[1], grid.shape[2], grid.shape[3]))
  else
    cartIndices = CartesianIndices(tuple(grid.shape...))
  end

  return [idx for idx in cartIndices if norm(collect(Tuple(idx)) .- grid.shape./2 .- 0.5) <= grid.shape[grid.radiusAxis]/2]
end

function getindex(grid::TubularRegularGridPositions, i::Integer)
  filteredPositions_ = filteredPositions(grid)
  if i > length(filteredPositions_) || i < 1
    throw(BoundsError(grid, i))
  end

  idx = Tuple(filteredPositions_[i])

  return ((-shape(grid).+(2 .*idx.-1))./shape(grid)).*fieldOfView(grid)./2 + fieldOfViewCenter(grid)
end

function getindex(grid::TubularRegularGridPositions, idx::Vector{T}) where T<:Number
  filteredPositions_ = filteredPositions(grid)
  cartIdx = CartesianIndex(idx...)
  if !(cartIdx in filteredPositions_)
    throw(BoundsError(grid, idx))
  else
    return ((-shape(grid).+(2 .*idx.-1))./shape(grid)).*fieldOfView(grid)./2 + fieldOfViewCenter(grid)
  end
end

function posToIdxFloat(grid::TubularRegularGridPositions, pos)
  idx = 0.5 .* (shape(grid) .* ((pos .- fieldOfViewCenter(grid)) ./
              ( 0.5 .* fieldOfView(grid) ) .+ 1) .+ 1)
  idx = [isnan(val) ? one(eltype(idx)) : val for val in idx]
  return idx
end

function posToIdx(grid::TubularRegularGridPositions, pos)
  return round.(Int64, posToIdxFloat(grid, pos))
end

function posToLinIdx(grid::TubularRegularGridPositions, pos)
  filteredPositions_ = filteredPositions(grid)
  idx = posToIdx(grid, pos)
  matchingIdx = [i for (i, idx_) in enumerate(filteredPositions_) if idx_ == CartesianIndex(idx...)]
  if length(matchingIdx) < 1
    error("The position $pos is not valid.")
  else
    return first(matchingIdx)
  end
end


# General functions for handling grids
fieldOfView(grid::GridPositions) = grid.fov
fieldOfView(grid::UniformRandomPositions{T, D, <:AxisAlignedBox}) where {T, D} = grid.domain.fov
fieldOfView(mgrid::MeanderingGridPositions) = fieldOfView(mgrid.grid)
fieldOfView(bgrid::BreakpointPositions) = fieldOfView(bgrid.grid)
shape(grid::GridPositions) = grid.shape
shape(mgrid::MeanderingGridPositions) = shape(mgrid.grid)
shape(bgrid::BreakpointPositions) = shape(bgrid.grid)
fieldOfViewCenter(grid::GridPositions) = grid.center
fieldOfViewCenter(grid::UniformRandomPositions) = grid.domain.center
fieldOfViewCenter(mgrid::MeanderingGridPositions) = fieldOfViewCenter(mgrid.grid)
fieldOfViewCenter(bgrid::BreakpointPositions) = fieldOfViewCenter(bgrid.grid)

spacing(grid::GridPositions) = grid.fov ./ grid.shape

struct SphericalTDesign{S, D, N, EL} <: Positions{S, D}
  T::UInt64
  radius::S
  positions::SMatrix{D, N, EL}
  center::SVector{D, S}
end

function SphericalTDesign(file::HDF5.File)
  T = read(file, "/positionsTDesignT")
  N = read(file, "/positionsTDesignN")
  radius = read(file, "/positionsTDesignRadius")*Unitful.m
  center = read(file, "/positionsCenter")*Unitful.m
  return loadTDesign(Int64(T),N,radius,center)
end
function SphericalTDesign(params::Dict)
  T = params["T"]
  N = params["N"]
  radius = params["radius"]*Unitful.m
  center = params["center"]*Unitful.m
  return loadTDesign(T, N, radius, center)
end

function write(file::HDF5.File, positions::SphericalTDesign)
  write(file,"/positionsType", "SphericalTDesign")
  write(file, "/positionsTDesignT", positions.T)
  write(file, "/positionsTDesignN", size(positions.positions,2))
  write(file, "/positionsTDesignRadius", Float64.(ustrip.(uconvert.(Unitful.m, positions.radius))) )
  write(file, "/positionsCenter", Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center)))) )
end
					
function toDict(positions::SphericalTDesign)
  params = Dict{String,Any}()
  params["type"] = "SphericalTDesign"
  params["T"] = positions.T
  params["N"] = size(positions.positions,2)
  params["radius"] = Float64.(ustrip.(uconvert.(Unitful.m, positions.radius)))
  params["center"] = Float64.(ustrip.(uconvert.(Unitful.m, Array(positions.center))))
  return params
end

getindex(tdes::SphericalTDesign, i::Integer) = tdes.radius.*tdes.positions[:,i] + tdes.center

const DEFAULT_TDESIGNS = @path joinpath(@__DIR__, "TDesigns.hd5")
"""
    loadTDesign(t::Int64, N::Int64, radius::S=10Unitful.mm, center::Vector{V}=[0.0,0.0,0.0]Unitful.mm, filename::String=joinpath(@__DIR__, "TDesigns.hd5")) where {S,V<:Unitful.Length}
*Description:* Returns the t-design array for chosen degree t and number of points N\\
\\
*Input:*
- `t` - degree
- `N` - number of points
- `radius` - radius of the sphere (default: 10.0mm)
- `center` - center of the sphere (default: [0.0,0.0,0.0]mm)
- `filename` - name of the file containing the t-designs (default loads TDesign.hd5)

*Output:*
- t-design of type SphericalTDesign in Cartesian coordinates containing t, radius, center and positions (which are located on the unit sphere unless `getindex(tdes,i)` is used)
"""
function loadTDesign(t, N, radius::S=10.00Unitful.mm, center::Vector{S}=[0.0,0.0,0.0]Unitful.mm, filename = DEFAULT_TDESIGNS) where {S<:Unitful.Length}
  h5file = h5open(filename, "r")
  address = "/$t-Design/$N"

  if haskey(h5file, address)
    positions = copy(transpose(read(h5file, address)))
    return SphericalTDesign(UInt(t),radius, SMatrix{3, size(positions, 2)}(positions), SVector{length(center)}(center))
  else
    if haskey(h5file, "/$t-Design/")
      Ns = Int[]
      for N in keys(read(h5file, string("/$t-Design")))
	push!(Ns,parse(Int,N))
      end
      sort!(Ns)
      @info "No spherical $t-Design with $N points available!\nThere are spherical $t-Designs with following N:" Ns
      throw(DomainError(1))
    else
      ts = Int[]
      for d in keys(read(h5file))
	m = match(r"(\d{1,})-(Design)",d)
	if m != nothing
	  push!(ts,parse(Int,m[1]))
        end
      end
      sort!(ts)
      @info "No spherical $t-Design available!\n Choose another t."
      throw(DomainError(1))
    end
  end
end

# Unstructured collection of positions
struct ArbitraryPositions{T, D, N} <: Positions{T, D}
  positions::SMatrix{D, N, T}
end
ArbitraryPositions(pos::AbstractMatrix) = ArbitraryPositions(SMatrix{size(pos, 1), size(pos, 2)}(pos)) 
getindex(apos::ArbitraryPositions, i::Integer) = apos.positions[:,i]

function ArbitraryPositions(grid::GridPositions{T, D}) where {T, D}
  positions = zeros(T,D,length(grid))
  for i=1:length(grid)
    positions[:,i] = grid[i]
  end
  return ArbitraryPositions(positions)
end

function write(file::HDF5.File, apos::ArbitraryPositions,)
  write(file,"/positionsType", "ArbitraryPositions")
  write(file, "/positionsPositions", Float64.(ustrip.(uconvert.(Unitful.m, Array(apos.positions)))) )
end

function ArbitraryPositions(file::HDF5.File)
  pos = read(file, "/positionsPositions")*Unitful.m
  return ArbitraryPositions(pos)
end

"""
    SortedPositions{T, D, G} <: NestedPositions{T, D, G}
    SortedPositions(grid::G, start::AbstractVector{T}=first(grid)) where {T, D, G <: Positions{T, D}}

Positions container which returns all points of the parent `grid` in a greedy nearest-neighbor order.
"""
struct SortedPositions{T, D, G} <: NestedPositions{T, D, G}
  parent::G
  indices::Vector{Int64}
  function SortedPositions(grid::G, start::AbstractVector{T} = first(grid)) where {T, D, G <: Positions{T, D}}
    current = start
    positions = collect(grid)
    sortedpos = Vector{typeof(start)}()
    indices = Int64[]
    # Greedily pick next position with smallest value
    while length(sortedpos) != length(grid)
      (val, idx) = findmin(map(x-> norm(x - current), positions))
      current = positions[idx]
      push!(sortedpos, current)
      push!(indices, idx)
      positions[idx] = [typemax(T) for i = 1:D]
    end
    return new{T, D, G}(grid, indices)
  end
end
length(grid::SortedPositions) = length(grid.indices)
parent(grid::SortedPositions) = grid.parent
parentindices(grid::SortedPositions) = grid.indices
getindex(grid::SortedPositions, i) = parent(grid)[grid.indices[i]]
fieldOfView(grid::SortedPositions) = fieldOfView(parent(grid))
fieldOfViewCenter(grid::SortedPositions) = fieldOfViewCenter(parent(grid))
shape(grid::SortedPositions) = shape(parent(grid))
spacing(grid::SortedPositions) = spacing(parent(grid))

# TODO: Specialize to make it faster
function Base.:(==)(val1::Positions, val2::Positions)
  for (x, y) in zip(val1, val2)
    if !(upreferred.(x) == upreferred.(y))
      return false
    end
  end

  return true
end

# fuction related to looping
length(tdes::SphericalTDesign) = size(tdes.positions,2)
length(apos::ArbitraryPositions) = size(apos.positions,2)
length(grid::GridPositions) = prod(grid.shape)
length(rpos::UniformRandomPositions) = rpos.N
length(mgrid::MeanderingGridPositions) = length(mgrid.grid)
length(bgrid::BreakpointPositions) = length(bgrid.grid)+length(bgrid.breakpointIndices)

start_(grid::Positions) = 1
next_(grid::Positions,state) = (grid[state],state+1)
done_(grid::Positions,state) = state > length(grid)
iterate(grid::Positions, s=start_(grid)) = done_(grid, s) ? nothing : next_(grid, s)
eltype(::Positions{T, D}) where {T, D} = SVector{D, T}

include("Interpolation.jl")
include("Subsampling.jl")