# MDF reconstruction data is loaded/stored using ImageMeta objects from
# the ImageMetadata.jl package.

using Images


export imcenter

imcenter(img::AxisArray) = map(x->(0.5*(last(x)+first(x))), ImageAxes.filter_space_axes(Images.axes(img), axisvalues(img)))
imcenter(img::ImageMeta) = imcenter(data(img))


    # # TODO: move the following to Analyze???
    # dateStr, timeStr = split("$(acqDate(b))","T")
    # dateStr = prod(split(dateStr,"-"))
    # timeStr = split(timeStr,".")[1]
    # timeStr = prod(split(timeStr,":"))
    #
    # header["date"] = dateStr
    # header["time"] = timeStr


function saveRecoDataMDF(filename, image::ImageMeta)
  L = size(image,ndims(image))

  C = size(image,1)
  N = div(length(data(image)), L*C)
  c = reshape(convert(Array,image), C, N, L )
  grid = size(image)[2:4]

  params = properties(image)
  params["recoData"] = c
  params["recoFov"] = collect(grid) .* collect(pixelspacing(image))
  params["recoFovCenter"] = collect(imcenter(image))
  params["recoSize"] = collect(grid)
  params["recoOrder"] = "xyz"

  h5open(filename, "w") do file
    saveasMDF(file, params)
  end

  if haskey(image,"recoParams")
    saveParams(file, "/reconstruction/parameters", image["recoParams"])
  end
end

function loadRecoDataMDF(filename::AbstractString)

  f = MDFFile(filename)

  header = loadMetadata(f)

  header["datatype"] = "MPI"
  pixspacing = recoFov(f) ./ recoSize(f)

  h5open(filename, "r") do file
    g = file["/reconstruction"]
    if exists(g, "parameters")
      params = loadParams(filename, "/reconstruction/parameters")
      header["recoParams"] = params
    end
  end

  c_ = h5read(f.filename, "/reconstruction/data")

  c = reshape(c_, size(c_,1), recoGridSize(f)..., size(c_,3))

  off = fieldOfViewCenter(f)
  if off != nothing
    offset = off .- 0.5.*recoFov(f) .+ 0.5.*pixspacing
  else
    offset = [0.0,0.0,0.0]

  im = AxisArray(c, (:color,:x,:y,:z,:time),
                      tuple(1.0, pixspacing..., dfcycle(f)),
                      tuple(0.0, offset..., 0.0))

  imMeta = ImageMeta(im, header)

  return imMeta
end

function saveParams(filename::AbstractString, path, params::Dict)
  h5open(filename, "w") do file
    saveParams(file, path, params)
  end
end

type ColoringParams
  cmin
  cmax
  cmap
end

function saveParams(file, path, params::Dict)
  for (key,value) in params
    ppath = joinpath(path,string(key))
    if typeof(value) <: Bool
      write(file, ppath, UInt8(value))
      dset = file[ppath]
      attrs(dset)["isbool"] = "true"
    elseif typeof(value) <: Range
      write(file, ppath, [first(value),step(value),last(value)])
      dset = file[ppath]
      attrs(dset)["isrange"] = "true"
    elseif value == nothing
      write(file, ppath, "")
      dset = file[ppath]
      attrs(dset)["isnothing"] = "true"
    elseif typeof(value) <: ColoringParams
      tmp = zeros(3)
      tmp[1] = value.cmin
      tmp[2] = value.cmax
      tmp[3] = value.cmap

      write(file, ppath, tmp)
      dset = file[ppath]
      attrs(dset)["iscoloring"] = "true"
    elseif typeof(value) <: Array{ColoringParams,1}
      tmp = zeros(3,length(value))
      for i=1:length(value)
        tmp[1,i] = value[i].cmin
        tmp[2,i] = value[i].cmax
        tmp[3,i] = value[i].cmap
      end
      write(file, ppath, tmp)
      dset = file[ppath]
      attrs(dset)["iscoloringarray"] = "true"
    elseif typeof(value) <: Array{Any}
      write(file, ppath, [v for v in value])
    else
      write(file, ppath, value)
    end
  end
end

function loadParams(filename::AbstractString, path)
  fid = h5open(filename, "r")
  params = loadParams(fid, path)
  close(fid)
  return params
end

function loadParams(file, path)
  params = Dict{Symbol,Any}()

  g = file[path]
  for obj in g
    key = last(splitdir(HDF5.name(obj)))
    data = read(obj)
    attr = attrs(obj)
    if exists(attr, "isbool")
      params[Symbol(key)] = Bool(data)
    elseif exists(attr, "isrange")
      if data[2] == 1
        params[Symbol(key)] = data[1]:data[3]
      else
        params[Symbol(key)] = data[1]:data[2]:data[3]
      end
    elseif exists(attr, "isnothing")
       params[Symbol(key)] = nothing
    elseif exists(attr, "iscoloring")
       params[Symbol(key)] = ColoringParams(data[1], data[2],round(Int64,data[3]))
    elseif exists(attr, "iscoloringarray")
       coloring = ColoringParams[]
       for i=1:size(data,2)
         push!(coloring, ColoringParams(data[1,i], data[2,i],round(Int64,data[3,i])))
       end
       params[Symbol(key)] = coloring
    else
      params[Symbol(key)] = data
    end
  end

  return params
end
