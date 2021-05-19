export changeParam

# Support for handling complex datatypes in HDF5 files
function writeComplexArray(file, dataset, A::AbstractArray{Complex{T},D}) where {T,D}
  d_type_compound = HDF5.h5t_create(HDF5.H5T_COMPOUND,2*sizeof(T))
  HDF5.h5t_insert(d_type_compound, "r", 0 , HDF5.hdf5_type_id(T))
  HDF5.h5t_insert(d_type_compound, "i", sizeof(T) , HDF5.hdf5_type_id(T))

  shape = collect(reverse(size(A)))
  space = HDF5.h5s_create_simple(D, shape, shape)

  dset_compound = HDF5.h5d_create(file, dataset, d_type_compound, space,
                                  HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT)
  HDF5.h5s_close(space)

  HDF5.h5d_write(dset_compound, d_type_compound, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, A)

  HDF5.h5d_close(dset_compound)
  HDF5.h5t_close(d_type_compound)
end

function isComplexArray(file, dataset)
  if eltype(file[dataset]) <: Complex
    return true
  else
    return false
  end
end

function getComplexType(file, dataset)
  T = HDF5.get_jl_type(
            HDF5.Datatype(
              HDF5.h5t_get_member_type( datatype(file[dataset]).id, 0 )
          )
        )
    return Complex{T}
end

function readComplexArray(file::HDF5.File, dataset)
  T = getComplexType(file, dataset)
  A = copy(HDF5.readmmap(file[dataset],getComplexType(file,dataset)))
  return A
end

function readComplexArray(filename::AbstractString, dataset)
  h5open(filename, "r") do file
    return readComplexArray(file, dataset)
  end
end

function changeParam(filename::AbstractString, paramName::AbstractString, paramValue)
  #GC.gc()
  h5open(filename, "r+") do file
	  if haskey(file, paramName)
      delete_object(file, paramName)
    end
    write(file, paramName, paramValue)
  end
end

function makeAxisArray(array::Array{T,5}, pixelspacing, offset, dt) where T
  N = size(array)
  im = AxisArray(array, Axis{:color}(1:N[1]),
		 Axis{:x}(range(offset[1],step=pixelspacing[1],length=N[2])),
		 Axis{:y}(range(offset[2],step=pixelspacing[2],length=N[3])),
		 Axis{:z}(range(offset[3],step=pixelspacing[3],length=N[4])),
		 Axis{:time}(range(0*unit(dt),step=dt,length=N[5])))
  return im
end

# RawFile removed see commit https://github.com/MagneticParticleImaging/MPIFiles.jl/commit/c2b49833cc00127770fb2e1f5342883f213ff002

"Macro for making required fields return `missing` instead of a key error."
macro keyrequired(expr)
  return quote
    try
      $(esc(expr))
    catch e
      if isa(e, KeyError)
        missing
      else
        rethrow()
      end
    end
  end
end

"Macro for making optional fields return `nothing` instead of a key error."
macro keyoptional(expr)
  return quote
    try
      $(esc(expr))
    catch e
      if isa(e, KeyError)
        nothing
      else
        rethrow()
      end
    end
  end
end
