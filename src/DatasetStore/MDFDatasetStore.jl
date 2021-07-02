export MDFDatasetStore, MDFStore, addStudy, getMDFStudyFolderName, calibdir, getMDFStore

struct MDFDatasetStore <: DatasetStore
  path::String

  function MDFDatasetStore(path::String)
    path = abspath(path)
    mkpath(path)
    mkpath(joinpath(path,"measurements"))
    mkpath(joinpath(path,"reconstructions"))
    mkpath(joinpath(path,"calibrations"))
    return new(path)
  end
end

if ispath("/opt/mpidata/Bruker")
  const MDFStore = MDFDatasetStore("/opt/data/Bruker")
end
path(e::Experiment{MDFDatasetStore}) = joinpath( path(e.study), string(e.num)*".mdf" )
path(s::Study{MDFDatasetStore}, numExp::Integer) = joinpath(path(s),string(numExp)*".mdf")
iscalib(e::Experiment{MDFDatasetStore}) = e.study.foldername == 
                             ".."*Base.Filesystem.path_separator*"calibrations"
readonly(::MDFDatasetStore) = false
changeParam(e::Experiment{MDFDatasetStore}, paramName::AbstractString, paramValue) =
               changeParam(path(e), paramName, paramValue)

studydir(d::MDFDatasetStore) = joinpath(d.path,"measurements")
calibdir(d::MDFDatasetStore) = joinpath(d.path,"calibrations")

function getMDFStore(study::Study{MDFDatasetStore})
    return study.store
end

function getMDFStore(experiment::Experiment{MDFDatasetStore})
    return getMDFStore(experiment.study)
end

function getStudy(d::MDFDatasetStore, studyfolder::String)
  study = nothing
  studypath = joinpath( studydir(d), studyfolder)
  if length(studyfolder) >= 15 &&
     isascii(studyfolder[1:15]) &&
     all([tryparse(Int,studyfolder[l:l])!=nothing for l=union(1:8,10:15)])

    w = split(studyfolder,'_')
    dateStr = w[1]
    timeStr = w[2]
    date = DateTime(string(dateStr[1:4],"-",dateStr[5:6],"-",dateStr[7:8],"T",
			   timeStr[1:2],":",timeStr[3:4],":",timeStr[5:6]))
    name = join(w[3:end],"_")
  else
    date = Dates.unix2datetime(stat(studypath).mtime)
    name = studyfolder
  end

  subject = ""
  study = Study(d, name, studyfolder, date, subject )
  return study
end

function getCalibStudy(d::MDFDatasetStore)
  return Study(d, "calibrations"; 
       foldername=".."*Base.Filesystem.path_separator*"calibrations") 
end

function getExperiments(s::Study{MDFDatasetStore})
  files = readdir(path(s))

  experiments = Experiment[]

  @debug "Time for get Experiments"
  for file in files
    prefix, ext = splitext(file)
    if !isdir(file) && tryparse(Int64,prefix) != nothing &&
       (ext == ".mdf" || ext == ".hdf" || ext == ".h5") &&
       isfile(joinpath(path(s),file))
       num = tryparse(Int64, prefix)
       if num != nothing 
        exp = getExperiment(s, num)
        push!(experiments, exp)
       end
    end
  end
  sort!(experiments,lt=(a,b)->(a.num < b.num))
  return experiments
end

getMDFStudyFolderName(study::Study) = getMDFStudyFolderName(study.name, study.date)

function getMDFStudyFolderName(name::String, date::DateTime)
  return string(split(string(date),"T")[1][union(1:4,6:7,9:10)],"_",
                split(string(date),"T")[2][union(1:2,4:5,7:8)],"_",name)
end

function addStudy(d::MDFDatasetStore, study::Study)
  studypath = joinpath( studydir(d), getMDFStudyFolderName(study))
  mkpath(studypath)
  try_chmod(studypath, 0o770, recursive=true)

  nothing
end

function Base.empty!(d::MDFDatasetStore)
  studies = getStudies(d)
  remove.(studies)
  rm(calibdir(d), recursive=true)
  mkpath(calibdir(d))
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
