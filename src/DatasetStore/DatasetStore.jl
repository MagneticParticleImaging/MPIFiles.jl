export Study, Experiment, Reconstruction, Visualization, DatasetStore,
       studydir, BrukerDatasetStore, BrukerStore, getStudy, getStudies, getExperiment,
       getExperiments, MDFDatasetStore, MDFStore, addReco, getReco, getRecons, findReco,
       findBrukerFiles, id, getVisus, getVisuPath, remove, addStudy, getNewExperimentNum,
       exportData, generateSFDatabase, loadSFDatabase, addVisu, readonly, getNewCalibNum,
       calibdir, try_chmod, getMDFStudyFolderName, path, makeTarGzip

include("Utils.jl")

abstract type DatasetStore end

# The following are base types describing
# the dataset store at a certain level
struct Study{D<:DatasetStore} 
  store::D
  name::String
  foldername::String
  date::DateTime
  subject::String
end

function Study(store::DatasetStore, name::String; foldername::String="", 
         date::DateTime = now(), subject::String="") 
  #remove Milliseconds
  date_ = DateTime(split(string(date),".")[1])
  if foldername == ""
    foldername = getMDFStudyFolderName(name,date_)
  end
  path = joinpath( studydir(store), foldername )
  mkpath(path)
  return Study(store, name, foldername, date_, subject)
end

path(s::Study) = joinpath( studydir(s.store), s.foldername )
id(s::Study) = s.name

# might not be so clever to use an explicit type fields here
# maybe better a dict
struct Experiment{D<:DatasetStore}
  study::Study{D}
  num::Int64
  name::String
  numFrames::Int64
  df::Vector{Float64}
  sfGradient::Float64
  numAverages::Int64
  operator::String
  time::String
  # more ...
end

mutable struct Reconstruction
  path::String
  num::Int64
  params::Dict
end

mutable struct Visualization
  path::String
  num::Int64
  params::Dict
end

include("BrukerDatasetStore.jl")
include("MDFDatasetStore.jl")
include("SFDatabase.jl")



### generic functions ###

function getStudies(d::DatasetStore)
  s = Study[]

  files = readdir( studydir(d) )
  for file in files
    fullpath = joinpath(studydir(d),file)
    if isdir(fullpath) && !ishidden(fullpath)
      try
      study = getStudy(d, file)
      if study != nothing
        push!(s, study)
      end
      catch e
        @warn e
      end
    end
  end
  return s
end

function remove(study::Study)
  if isdir(path(study))
    rm(path(study), recursive=true)

    #TODO remove recos!
  end
end

function getExperiment(s::Study, numExp::Integer)

  path_ = path(s, numExp)

  if ispath(path_)
    b = MPIFile(path_, fastMode=true)
  else
    return nothing
  end
  prefix, ext = splitext(path_)

  exp = Experiment(s, parse(Int64,last(splitdir(prefix))), #why parse experiment number here????, 
                      string(experimentName(b)), acqNumFrames(b),
                      round.(1000 .* vec(dfStrength(b)[1,:,1]),digits=2), maximum(abs.(acqGradient(b))),
                      acqNumAverages(b), scannerOperator(b), string(acqStartTime(b)))

  return exp
end

function remove(exp::Experiment)
  if isfile(exp.path)
    rm(exp.path)
  end
end

function exportData(e::Experiment, mdf::MDFDatasetStore; logging=false, storeForwardRef=false, kargs...)
  # pretend to be a measurement to enforce loading data from time domain in case post processed data is not available
  # in case of compression isCalib=false right now does not work, thus we disable it.
  
  isCalib = (haskey(kargs, :SNRThresh) || haskey(kargs, :sparsityTrafoRedFactor)) ? true : false

  f = MPIFile(path(e), isCalib=isCalib)

  if !isConvertibleToMDF(f)
    return ""
  end

  if iscalib(e)
    exportpath = getNewCalibPath(mdf)
    saveasMDF(exportpath, f; applyCalibPostprocessing=true, kargs...)
    @info "Calibration data from $path successfully exported to $exportpath." 
  else
    name = studyName(f)
    subject = experimentSubject(f)
    date = studyTime(f)
    s = Study(mdf, name; date=date, subject=subject)
    exportpath = getNewExperimentPath(s)
    saveasMDF(exportpath, f; kargs...)
    @info "Measurement data from $path successfully exported to $exportpath." 
  end

  if isdir(path(e)) && storeForwardRef # only if input data is a BrukerFile
    # Store export path in Bruker directory
    open(joinpath(path(e),"mdf"),write=true) do io
      write(io, exportpath)
    end
  end

  if logging
    # log action
    open("/opt/DataArchiveOptmpidata/convertBrukerToMDF/log.csv", append=true) do io
      write(io,"$path, $exportpath\n")
    end
  end

  return exportpath
end

function exportData(s::Study, mdf::MDFDatasetStore; kargs...)
  exps = getExperiments(s)
  for e in exps
    exportData(e, mdf; kargs...)
  end
end

function exportData(store::D, mdf::MDFDatasetStore; kargs...) where D<:DatasetStore
  for s in getStudies(store)
    addStudy(mdf, s)
    exportData(s, mdf; kargs...)
  end

  if D == MDFDatasetStore
    exportData(getCalibStudy(store), mdf; kargs...) 
  end
end

#########

# generate tar balls and Julia Artifacts.
function makeTarGzip(dir, filename)
  tar_gz = open(filename, write=true)
  tar = GzipCompressorStream(tar_gz)
  Tar.create(dir, tar)
  close(tar)
end

function makeTarGzip(d::DatasetStore, filename)
  makeTarGzip(d.path, filename)
end

function makeTarGzip(d::DatasetStore)
  filename = joinpath(splitpath(d.path)...)*".tar.gz" # this ensures that a trailing / is removed
  makeTarGzip(d.path, filename)
end



#########


function getNewNumInFolder(path)
  if !isdir(path)
    mkpath(path)
    try_chmod(path, 0o777, recursive=true)
    return 1
  end

  files = readdir(path)
  num = 1
  if length(files) > 0
    for i=1:length(files)
      pref, ext = splitext(files[i])
      num_ = tryparse(Int64, pref)
      if num_ != nothing && num_+1>num
        num = num_+1
      end
    end
  end

  return num
end

function getNewExperimentNum(s::Study)
  return getNewNumInFolder(path(s))
end

function getNewExperimentPath(s::Study)
  addStudy(s.store, s)
  expNum = getNewExperimentNum(s)
  path_ = joinpath(path(s),string(expNum)*".mdf")
  # touch new mdf file
  touch(path_)
  try_chmod(path_, 0o660)
  return path_
end

function getNewCalibNum(d::MDFDatasetStore)
  return getNewNumInFolder(calibdir(d))
end

function getNewCalibPath(d::MDFDatasetStore)
    calibNum = getNewCalibNum(d)
    path = joinpath(calibdir(d),string(calibNum)*".mdf")
    # touch new mdf file
    touch(path)
    try_chmod(path, 0o660)
    return path
end

include("Reconstruction.jl")
include("Visualization.jl")
