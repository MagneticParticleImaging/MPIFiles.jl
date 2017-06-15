import FileIO: save

export Study, Experiment, Reconstruction, Visualization, DatasetStore,
       studydir, BrukerDatasetStore, BrukerStore, getStudy, getStudies, getExperiment,
       getExperiments, MDFDatasetStore, MDFStore, addReco, getReco, getRecons, findReco,
       findBrukerFiles, id, getVisus, getVisuPath, remove

########################################

# The following are base types describing
# the dataset store at a certain level

type Study
  path::String
  name::String
  subject::String
  date
end

id(s::Study) = s.name

# might not be so clever to use explicity type fields here
# maybe better a dict
type Experiment
  path::String
  num::Int64
  description::String
  numFrames::Int64
  dfFov::Vector
  sfGradient::Vector
  numAverages::Int64
  operator::String
  # more ...
end

type Reconstruction
  path::String
  num::Int64
  params::Dict
end

type Visualization
  path::String
  num::Int64
  params::Dict
end

abstract DatasetStore

type BrukerDatasetStore <: DatasetStore
   path::String
end

const BrukerStore = BrukerDatasetStore("/opt/mpidata")

type MDFDatasetStore <: DatasetStore
  path::String
end

const MDFStore = MDFDatasetStore("/opt/data")

### generic functions ###

function ishidden(filename::AbstractString)
  @static if is_unix()
    s = basename(filename)
    return (!isempty(s) && s[1] == '.')
  else
    attr = ccall((:GetFileAttributesA), stdcall, Cint, (Ptr{Uint8},),bytestring(filename))
    return attr & 0x2 > 0
  end
end

function getStudies(d::DatasetStore)
  s = Study[]

  files = readdir( studydir(d) )
  for file in files
    fullpath = joinpath(studydir(d),file)
    if isdir(fullpath) && !ishidden(fullpath)
      study = getStudy(d, file)
      if study != nothing
        push!(s, study)
      end
    end
  end
  return s
end

function getExperiment(path::String)

  prefix, ext = splitext(path)

  if isdir(path) #Ugly
    b = MPIFiles.BrukerFileFast(path) #use fast path for BrukerFiles
  else
    b = MDFFile(string(prefix,".mdf"))
  end

  exp = Experiment( path, parse(Int64,last(splitdir(prefix))),
                      string(experimentDescription(b)), measNumFrames(b),
                      round(1000.*vec(acqFov(b)),2), acqGradient(b)[:,1],
                      rxNumAverages(b), scannerOperator(b))

  return exp
end

getExperiment(s::Study, numExp::Integer) = getExperiment(joinpath(s.path,string(numExp)))

function exportToMDFStore(d::BrukerDatasetStore, s::Study, e::Experiment, mdf::MDFDatasetStore, freqSpace=false)

  name = s.name*"_MDF"
  path = joinpath( studydir(mdf), name)
  subject = s.subject
  date = ""

  newStudy = Study(path,name,subject,date)

  addStudy(mdf, newStudy)
  expNum = getNewExperimentNum(mdf, newStudy)

  b = BrukerFile(e.path)

  if isSF(b)
    calibNum = getNewCalibNum(mdf)
    #saveasMDF( joinpath(calibdir(mdf),string(calibNum)*".mdf"), b, deltaSampleSize=[0.002,0.002,0.001]) #Fixme
    saveasMDF( joinpath(calibdir(mdf),string(calibNum)*".mdf"),
               b, deltaSampleSize=[0.0,0.0,0.0], bgcorrection=true ) #Fixme
  else
    saveasMDF( joinpath(studydir(mdf),newStudy.name,string(expNum)*".mdf"), b, timeDomain=!freqSpace)
  end

end

###  Implementations of abstract interfaces ###

studydir(d::BrukerDatasetStore) = d.path
studydir(d::MDFDatasetStore) = joinpath(d.path,"measurements")
calibdir(d::MDFDatasetStore) = joinpath(d.path,"calibrations")

function getStudy(d::BrukerDatasetStore, studyfolder::String)
  study = nothing
  studypath = joinpath(d.path,studyfolder)
  if !ishidden(studypath) && isdir(studypath)
    w = split(studyfolder,'_')
    if length(w) >= 5 # only these can be study folders
      w_ = w[1:end-2]
      date = w[1]
      date = string(date[1:4],"/",date[5:6],"/",date[7:8])


      j = JcampdxFile()
      subjfile = string(studypath,"/subject")
      if isfile(subjfile)
        read(j,string(studypath,"/subject"),maxEntries=14) #magic number...
        name = string(latin1toutf8(j["SUBJECT_name_string"]),
                      "_",latin1toutf8(j["SUBJECT_study_name"]),
                      "_",latin1toutf8(j["SUBJECT_study_nr"]))
        subject = latin1toutf8(j["SUBJECT_name_string"])
      else
        # Workaround if no subject file is present => use first dataset
        # and derive the study from the Brukerfile
        r = readdir(studypath)

        found = false
        for file in r

          if !isnull(tryparse(Int64,file))

            b = BrukerFile(joinpath(studypath, file ))
            name = studyName(b)
            subject = experimentSubject(b)
            found = true
            break
          end
        end
        if !found
          return nothing
        end
      end
      study = Study(studypath, name, subject, date )
    end
  end
  return study
end

function getStudy(d::MDFDatasetStore, studyfolder::String)
  study = nothing
  studypath = joinpath( studydir(d), studyfolder)
  name = studyfolder
  subject = ""
  date = split(string(Dates.unix2datetime(stat(studypath).mtime)),"T")[1]
  study = Study(studypath, name, subject, date )
  return study
end

function addStudy(d::MDFDatasetStore, study::Study)
  studypath = joinpath( studydir(d), study.name)
  mkpath(studypath)

  nothing
end

function findBrukerFiles(path::AbstractString)
  files = readdir(path)

  bfiles = String[]

  for file in files
    if isdir(joinpath(path,file))
     try
      if isfile(joinpath(path,file,"acqp"))
        push!(bfiles, joinpath(path,file))
      else
        rfiles = findBrukerFiles(joinpath(path,file))
        if rfiles != nothing && length(rfiles) > 0
          push!(bfiles, rfiles...)
        end
      end
     catch
      continue
     end
    end
  end
  bfiles
end

function getExperiments(d::BrukerDatasetStore, s::Study)

  files = findBrukerFiles(s.path) # make me fast

  experiments = Experiment[]

  for file in files
    try
      exp = getExperiment(file)

      push!(experiments, exp)
    catch e
      println(e)
    end
  end
  return experiments
end

function getExperiments(d::MDFDatasetStore, s::Study)

  files = readdir(s.path)

  experiments = Experiment[]

  for file in files
    prefix, ext = splitext(file)
    if !isdir(file) && !isnull(tryparse(Int64,prefix)) &&
       (ext == ".mdf" || ext == ".hdf" || ext == ".h5")
      try
        exp = getExperiment(joinpath(s.path,file))

        push!(experiments, exp)
      catch e
        println(e)
      end
    end
  end
  return experiments
end

function getNewNumInFolder(d::MDFDatasetStore, path)
  if !isdir(path)
    mkpath(path)
    return 1
  end

  files = readdir(path)
  num = 1
  if length(files) > 0
    for i=1:length(files)
      pref, ext = splitext(files[i])
      num_ = tryparse(Int64, pref)
      if !isnull(num_) && get(num_)+1>num
        num = get(num_)+1
      end
    end
  end

  return num
end

function getNewExperimentNum(d::MDFDatasetStore, s::Study)
  return getNewNumInFolder(d, s.path)
end

function getNewCalibNum(d::MDFDatasetStore)

  files = readdir(calibdir(d))
  calibNum = 1
  if length(files) > 0
    for i=1:length(files)
      pref, ext = splitext(files[i])
      num = parse(Int64, pref) + 1
      if num>calibNum
        calibNum = num
      end
    end
  end

  return calibNum
end

####### Reconstruction Store MDF ###################

function getReco(d::MDFDatasetStore, study::Study, exp::Experiment, recoNum::Int64)
  path = joinpath(d.path, "reconstructions", id(study), string(exp.num), string(recoNum))
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

  outputpath = joinpath(d.path, "reconstructions", id(study), string(exp.num))
  # create data directory
  mkpath(outputpath)

  recoNum = getNewNumInFolder(d, outputpath)

  filepath = joinpath(outputpath, string(recoNum))

  saveRecoDataMDF(filepath*".mdf", image)
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
    if exists(file, "/reconstruction/parameters")
      o_delete(file, "/reconstruction/parameters")
    end
    saveParams(file, "/reconstruction/parameters", reco.params)
  end
end

function loadParams(reco::Reconstruction)

  if isfile(reco.path)
   h5open(reco.path, "r") do file
    g = file["/reconstruction"]
    if exists(g, "parameters") #new world order
      reco.params = loadParams(reco.path, "/reconstruction/parameters")
    else #this needs to go
      println("opening legacy file")
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

  datadir = joinpath(d.path, "reconstructions", id(study), string(exp.num))

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















####### Visualization Store #######


function getVisu(d::MDFDatasetStore, study::Study, exp::Experiment, reco::Reconstruction, numVisu)

  filename = joinpath(d.path, "reconstructions", id(study), string(exp.num), string(reco.num)*".visu")

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

  filename = joinpath(d.path, "reconstructions", id(study), string(exp.num), string(reco.num)*".visu")

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
    if exists(file, string(visu.num))
      o_delete(file, string(visu.num))
    end
  end
end

function save(visu::Visualization)
  if isfile(visu.path)
    file = h5open(visu.path, "r+")
  else
    file = h5open(visu.path, "w")
  end
  if exists(file, string(visu.num))
    o_delete(file, string(visu.num))
  end
  saveParams(file, string(visu.num), visu.params)
  close(file)
end

function addVisu(d::MDFDatasetStore, study::Study, exp::Experiment, reco::Reconstruction, visuParams)

  filename = joinpath(d.path, "reconstructions", id(study), string(exp.num), string(reco.num)*".visu" )

  visus = getVisus(d, study, exp, reco)

  if isempty(visus)
    num = 1
  else
    num = last(visus).num + 1
  end

  visu =  Visualization(filename, num, visuParams)

  save(visu)
end
