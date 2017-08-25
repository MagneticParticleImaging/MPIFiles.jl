# MDF reconstruction data is loaded/stored using ImageMeta objects from
# the ImageMetadata.jl package.

using Images


export imcenter, loadRecoDataMDF, saveRecoDataMDF

imcenter(img::AxisArray) = map(x->(0.5*(last(x)+first(x))), ImageAxes.filter_space_axes(Images.axes(img), axisvalues(img)))
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
  params["recoFovCenter"] = collect(imcenter(image))
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

  header = loadMetadata(f)

  header["datatype"] = "MPI"
  pixspacing = recoFov(f) ./ recoSize(f)

  recoParams = recoParameters(f)
  if recoParams != nothing
    header["recoParams"] = recoParams
  end

  c_ = recoData(f)

  c = reshape(c_, size(c_,1), recoSize(f)..., size(c_,3))

  off = recoFovCenter(f)
  if off != nothing
    offset = off .- 0.5.*recoFov(f) .+ 0.5.*pixspacing
  else
    offset = [0.0,0.0,0.0]
  end
  im = AxisArray(c, (:color,:x,:y,:z,:time),
                      tuple(1.0, pixspacing..., acqFramePeriod(f)),
                      tuple(0.0, offset..., 0.0))

  imMeta = ImageMeta(im, header)

  return imMeta
end
