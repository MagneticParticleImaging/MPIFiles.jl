export Study, getStudy, getStudies

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
@mustimplement path(s::Study{D}, numExp::Integer) where D<:DatasetStore

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