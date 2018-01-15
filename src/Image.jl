# MDF reconstruction data is loaded/stored using ImageMeta objects from
# the ImageMetadata.jl package.

export imcenter, loadRecoDataMDF, saveRecoDataMDF

imcenter(img::AxisArray) = map(x->(0.5*(last(x)+first(x))), ImageAxes.filter_space_axes(ImageAxes.axes(img), axisvalues(img)))
imcenter(img::ImageMeta) = imcenter(data(img))


function saveRecoDataMDF(filename, image::ImageMeta)
  C = colordim(image) == 0 ? 1 : size(image,colordim(image))
  L = timedim(image) == 0 ? 1 : size(image,timedim(image))

  if colordim(image) == 0
    grid = size(image)[1:3]
  else
    grid = size(image)[2:4]
  end
  N = div(length(data(image)), L*C)
  c = reshape(convert(Array,image), C, N, L )

  params = properties(image)
  params["recoData"] = c
  params["recoFov"] = collect(grid) .* collect(pixelspacing(image))[1:3]
  params["recoFovCenter"] = collect(imcenter(image))[1:3]
  params["recoSize"] = collect(grid)
  params["recoOrder"] = "xyz"
  if haskey(params,"recoParams")
    params["recoParameters"] = params["recoParams"]
  end

  h5open(filename, "w") do file
    saveasMDF(file, params)
  end
end

function loadRecoDataMDF(filename::AbstractString)
  f = MDFFile(filename)
  return loadRecoDataMDF(f)
end

function loadRecoDataMDF(f::MDFFile)
  header = loadMetadata(f)
  header["datatype"] = "MPI"

  recoParams = recoParameters(f)
  if recoParams != nothing
    header["recoParams"] = recoParams
  end
  im = loadRecoDataMDF_(f)
  imMeta = ImageMeta(im, header)

  return imMeta
end

function loadRecoDataMDF_(f::MDFFile)

  rsize::Vector{Int64} = recoSize(f)

  pixspacing::Vector{Float64} = recoFov(f) ./ rsize

  c_::Array{Float32,3} = recoData(f)
  #c_::Array{Float32,3} = h5read(filename, "/reconstruction/data")
  c::Array{Float32,5} = reshape(c_, size(c_,1), rsize[1], rsize[2], rsize[3], size(c_,3))

  off::Vector{Float64} = vec(recoFovCenter(f))
  offset::Vector{Float64} = [0.0,0.0,0.0]
  if off != nothing
    offset[:] = off .- 0.5.*recoFov(f) .+ 0.5.*pixspacing
  end
  periodTime = Float64(acqFramePeriod(f))

  ax1 = Axis{:color}(range(0.0, 1.0, size(c,1)))
  ax2 = Axis{:x}(range(offset[1], pixspacing[1], size(c,2)))
  ax3 = Axis{:y}(range(offset[2], pixspacing[2], size(c,3)))
  ax4 = Axis{:z}(range(offset[3], pixspacing[3], size(c,4)))
  ax5 = Axis{:time}(range(0.0, periodTime, size(c,5)))
  im = AxisArray(c,ax1,ax2,ax3,ax4,ax5)

  return im
end

#precompile(MPIFiles.loadRecoDataMDF_,(MDFFileV1,))
#precompile(MPIFiles.loadRecoDataMDF_,(MDFFileV2,))
