export AbstractSubsampledPositions, SubsampledPositions
"""
    AbstractSubsampledPositions{T, D, G} <: Positions{T, D}

Abstract supertype for positions that represent a subsample (by linear indices) of a parent
`Positions` object. The subsample is a view-like selector

Type parameters:
- T: coordinate element type
- D: spatial dimensionality
- G: parent `Positions` concrete type
"""
abstract type AbstractSubsampledPositions{T, D, G} <: Positions{T, D} end
"""
    parent(grid::AbstractSubsampledPositions)

Return the parent `Positions` object from which this subsampled view selects its elements.

Must be implemented by subtypes of `AbstractSubsampledPositions`.
"""
parent(grid::AbstractSubsampledPositions) = error("$(typeof(grid)) must implement `parent`")
"""
    parentindices(grid::AbstractSubsampledPositions)

Return the indices in the parent which correspond to the subsampling of `grid``.

Must be implemented by subtypes of `AbstractSubsampledPositions`.
"""
parentindices(grid::AbstractSubsampledPositions, i) = error("$(typeof(grid)) must implement `parentindices`")
getindex(grid::AbstractSubsampledPositions, i) = parent(grid)[parentindices(grid)[i]]
fieldOfView(grid::AbstractSubsampledPositions) = fieldOfView(parent(grid))
fieldOfViewCenter(grid::AbstractSubsampledPositions) = fieldOfViewCenter(parent(grid))
# I think one could argue if these are valid or not
shape(grid::AbstractSubsampledPositions) = shape(parent(grid))
spacing(grid::AbstractSubsampledPositions) = spacing(parent(grid))

"""
    SubsampledPositions{T, D, G <: Positions{T, D}} <: AbstractSubsampledPositions{T, D, G}

A concrete positions container that holds a reference to a parent `Positions` object and a vector
of linear indices selecting a subset of the parent’s elements.
"""
struct SubsampledPositions{T, D, G <: Positions{T, D}} <: AbstractSubsampledPositions{T, D, G}
  parent::G
  indices::Vector{Int64}
  function SubsampledPositions(grid::G, indices::Vector{Int64}) where {T, D, G <: Positions{T, D}}
    if !all(idx -> 1 <= idx <= length(grid), indices)
      throw(ArgumentError("Provided subsampling indices don't match the underlying parent positions."))
    end
    return new{T, D, G}(grid, indices)
  end
end
"""
    SubsampledPositions(grid, factor::Float64; seed = rand(UInt64))

Create a subsample of `grid` by selecting `round(length(grid) * factor)` unique positions
at random using a stable RNG seeded with `seed`. The selection is a uniform shuffle of
`1:length(grid)` truncated to the requested count.

Notes:
- The returned indices are unique and in random order.
- The same `seed` yields identical indices for the same `grid` length/indices.
"""
SubsampledPositions(grid, factor::Float64; seed = rand(UInt64)) = SubsampledPositions(grid, round(UInt64, length(grid) * factor); seed = seed)
"""
    SubsampledPositions(grid, numIndices::Integer; seed = rand(UInt64))

Create a subsample of `grid` by selecting `numIndices` unique positions uniformly at random
(using a stable RNG seeded with `seed`) from `1:length(grid)`. The indices are stored in random order.
"""
function SubsampledPositions(grid, numIndices::Integer; seed = rand(UInt64))
  if numIndices > length(grid)
    throw(ArgumentError("Requested more random `SubsampledPositions` indices ($numIndices) than exist in the parent grid ($(length(grid)))"))
  end
  rng = StableRNG(seed)
  indices = shuffle(rng, 1:length(grid))[1:numIndices]
  return SubsampledPositions(grid, indices)
end
"""
    SubsampledPositions(grid, indices_iterable)

Create a subsample of `grid` with the specified iterable of linear indices. The indices
are collected into vector in the provided order.
"""
SubsampledPositions(grid, other) = SubsampledPositions(grid, collect(other))
"""
    SubsampledPositions(grid::Positions{T}, positions::AbstractMatrix{T}) where T

Construct a `SubsampledPositions` by locating each column of `positions` (a D×N matrix of
coordinates with element type `T`) in the iterable `grid::Positions{T}`. Each column is matched
by exact equality against the coordinates returned by iterating `grid`, and the resulting
subsample preserves the column order and allows duplicates.

Notes:
- Matching is exact (no tolerance). Coordinate types must match (`T`) and values must be equal.
- Order and duplicate columns are preserved in the resulting subsample.
"""
function SubsampledPositions(grid::Positions{T, D}, positions::AbstractMatrix{T}) where {T, D}
  if size(positions, 1) != D
    throw(ArgumentError("Dimension of grid $D does not match dimension of positions $(size(positions, 1))"))
  end

  # This method uses exact equality atm, one could also implement an approximate version.
  # That would be more expensive (O(n^2)) though
  # Enforce known "key" type
  helper = Dict{SVector{D, T}, Int64}(SVector{D}(pos) => i for (i, pos) in enumerate(grid))
  indices = fill(0, size(positions, 2))

  for (i, pos) in enumerate(eachcol(positions))
    # Use same "key" type
    posV = SVector{D}(pos)
    indices[i] = get(helper, posV, 0)
  end

  # All positions should have been found
  noMatchingIndex = findall(idx -> idx == 0, indices)
  if !isempty(noMatchingIndex)
    throw(ArgumentError("Found no matching position in the grid for positions in columns: $noMatchingIndex"))
  end

  return SubsampledPositions(grid, indices)
end
"""
    SubsampledPositions(shape, fov, center::AbstractVector{T}, positions::AbstractMatrix{T})

Convenience constructor that builds a `RegularGridPositions(shape, fov, center)` and then
creates a `SubsampledPositions` by locating each column of `positions` in that grid.
"""
SubsampledPositions(shape, fov, center::AbstractVector{T}, positions::AbstractMatrix{T}) where T = SubsampledPositions(RegularGridPositions(shape, fov, center), positions) 

view(grid::Positions, inds...) = SubsampledPositions(grid, inds...)
length(grid::SubsampledPositions) = length(grid.indices)
parent(grid::SubsampledPositions) = grid.parent
parentindices(grid::SubsampledPositions) = grid.indices