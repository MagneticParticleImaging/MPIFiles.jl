export Reconstruction, addReco, getReco, getRecons, findReco

# TODO: make store a member and remove path
mutable struct Reconstruction
  path::String
  num::Int64
  params::Dict
end

function getReco(d::MDFDatasetStore, study::Study, exp::Experiment, recoNum::Int64)
  path = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study), string(exp.num), string(recoNum))
  filename = path*".mdf"
  if !isfile(filename)
    filename = path*".hdf"
  end
  r = Reconstruction(filename, recoNum, Dict())
  loadParams(r)
  return r
end

# This functions searches for recoparams and returns the corresponding recoNumber
# 0 indicates that the set of parameters was not found
function findReco(d::MDFDatasetStore, study::Study, exp::Experiment, recoParams::Dict)
  recoNum = 0

  recoParams_ = deepcopy(recoParams)
  # We do not care if the reconstruction has been done by a different
  # user. Therefore, we remove the :reconstructor field
  if haskey(recoParams_, :reconstructor)
    delete!(recoParams_, :reconstructor)
  end
  recons = getRecons(d, study, exp)
  for reco in recons
    if haskey(reco.params, :reconstructor)
      delete!(reco.params, :reconstructor)
    end
    if recoParams_ == reco.params
      recoNum = reco.num
    end
  end

  return recoNum
end

function getReco(d::MDFDatasetStore, study::Study, exp::Experiment, recoParams::Dict)
  getReco(d, study, exp, findReco(d,study,exp,recoParams) )
end


# The following function is certainly not ideal when considering a "getReco" scenario
function addReco(d::MDFDatasetStore, study::Study, exp::Experiment, image)

  outputpath = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study), string(exp.num))
  # create data directory
  mkpath(outputpath)
  try_chmod(outputpath, 0o777, recursive=true)

  recoNum = getNewNumInFolder(outputpath)

  filepath = joinpath(outputpath, string(recoNum))

  saveRecoData(filepath*".mdf", image)
  #save(filepath*".jld","recoParams",recoParams)
end


function remove(reco::Reconstruction)
  if isfile(reco.path)
    rm(reco.path)
  end
  visufile = getVisuPath(reco)
  if isfile(visufile)
    rm(visufile)
  end
end

function save(reco::Reconstruction)
  h5open(reco.path, "r+") do file
    if haskey(file, "/reconstruction/_parameters")
      delete_object(file, "/reconstruction/_parameters")
    end
    saveParams(file, "/reconstruction/_parameters", reco.params)
  end
end

function loadParams(reco::Reconstruction)

  if isfile(reco.path)
   h5open(reco.path, "r") do file
    g = file["/reconstruction"]
    if haskey(g, "_parameters") #new world order
      reco.params = loadParams(reco.path, "/reconstruction/_parameters")
    else #this needs to go
      @debug "opening legacy file"
      prefix, ext = splitext(reco.path)
      reco.params = load(prefix*".jld","recoParams")
    end
   end
  end
  nothing
end

function getVisuPath(reco::Reconstruction)
  prefix, ext = splitext(reco.path)
  return prefix*".visu"
end


function getRecons(d::MDFDatasetStore, study::Study, exp::Experiment)

  recons = Reconstruction[]

  datadir = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study), string(exp.num))

  if isdir(datadir)
    files = readdir(datadir)
    for file in files
        prefix, ext = splitext(file)
        fullfile = joinpath(datadir,file)


        if ext == ".hdf"  || ext == ".mdf"

        filename = joinpath(datadir,prefix*".mdf")
        if !isfile(filename)
            filename = joinpath(datadir,prefix*".hdf")
        end

        num = parse(Int64,prefix)
        r = Reconstruction(filename, num,  Dict())
        loadParams(r)

        push!(recons, r )
        end
    end
  end

  return recons
end
