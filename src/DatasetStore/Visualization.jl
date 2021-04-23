export Visualization, getVisus, getVisuPath, addVisu

mutable struct Visualization
  path::String
  num::Int64
  params::Dict
end

function getVisu(d::MDFDatasetStore, study::Study, exp::Experiment, reco::Reconstruction, numVisu)

  filename = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study),
			       string(exp.num), string(reco.num)*".visu")

  if isfile(filename)

  #  file = h5open(filename, "r")
 #   g = file[string(numVisu)]
    params = loadParams(filename, string(numVisu))

#    close(file)
  else
    params = Dict{Symbol,Any}()
  end
  return Visualization(filename, numVisu, params)
end

function getVisus(d::MDFDatasetStore, study::Study, exp::Experiment, reco::Reconstruction)

  visus = Visualization[]

  filename = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study), string(exp.num),
			      string(reco.num)*".visu")

  if isfile(filename)

    file = h5open(filename, "r")
    g = file["/"]
    for obj in g
      key = HDF5.name(obj)
      numVisu = parse(Int64, last(splitdir(key)) )
      params = loadParams(file, key) #this is the dictionarry

      push!(visus, Visualization(filename, numVisu, params))
    end
    close(file)
  end
  return visus
end

function remove(visu::Visualization)
  h5open(visu.path, "r+") do file
    if haskey(file, string(visu.num))
      delete_object(file, string(visu.num))
    end
  end
end

function save(visu::Visualization)
  if isfile(visu.path)
    file = h5open(visu.path, "r+")
  else
    file = h5open(visu.path, "w")
  end
  if haskey(file, string(visu.num))
    delete_object(file, string(visu.num))
  end
  saveParams(file, string(visu.num), visu.params)
  close(file)
end

function addVisu(d::MDFDatasetStore, study::Study, exp::Experiment, reco::Reconstruction, visuParams)

  filename = joinpath(d.path, "reconstructions", getMDFStudyFolderName(study), string(exp.num),
			      string(reco.num)*".visu" )

  visus = getVisus(d, study, exp, reco)

  if isempty(visus)
    num = 1
  else
    num = last(visus).num + 1
  end

  visu =  Visualization(filename, num, visuParams)

  save(visu)
end
