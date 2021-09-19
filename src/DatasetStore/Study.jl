export Study, getStudy, getStudies, validate, enforceStudy, changeStudy

# The following are base types describing
# the dataset store at a certain level
struct Study{D<:DatasetStore} 
  store::D
  name::String
  foldername::String
  date::DateTime
  subject::String
  uuid::UUID
end

function Study(store::DatasetStore, name::String, foldername::String, date::DateTime, subject::String)
  rng = StableRNG(hash(foldername)) 
  uuid = uuid4(rng)
  return Study(store, name, foldername, date, subject, uuid)
end

function Study(store::DatasetStore, name::String; foldername::String="", 
         date::DateTime = now(), subject::String="", createDir=true) 
  #remove Milliseconds
  date_ = trunc(date, Dates.Second)
  if foldername == ""
    foldername = getMDFStudyFolderName(name,date_)
  end
  path = joinpath( studydir(store), foldername )
  if createDir
    mkpath(path)
  end
  rng = StableRNG(hash(foldername)) 
  uuid = uuid4(rng)
  return Study(store, name, foldername, date_, subject)
end

path(s::Study) = joinpath( studydir(s.store), s.foldername )
@mustimplement path(s::Study{D}, numExp::Integer) where D<:DatasetStore

function getStudies(d::DatasetStore, name::String="")
  s = Study[]

  files = readdir( studydir(d) )
  for file in files
    fullpath = joinpath(studydir(d),file)
    if isdir(fullpath) && !ishidden(fullpath)
      try
      study = getStudy(d, file)
      # return only studies where name is a substring of the studyname
      if study !== nothing && match(Regex("(?i)"*name), study.name) !== nothing
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


function changeParam(s::Study, paramName::AbstractString, paramValue) 
  exps = getExperiments(s)
  for e in exps
    changeParam(e, paramName, paramValue)
  end
end

function changeStudy(s::Study, newStudyName::AbstractString; date::DateTime = now()) 
  changeParam(s, "/study/name", newStudyName) 
  changeParam(s, "/study/time", string(date)) 
  sNew = Study(s.store, newStudyName, date=date, subject=s.subject, createDir=false)
  mv(path(s), path(sNew))
  return 
end

function validate(s::Study)
  exps = getExperiments(s)
  valid = true
  for e in exps
    MPIFile(path(e), fastMode=true) do f
      date1 = trunc(studyTime(f), Dates.Second)
      date2 = trunc(s.date, Dates.Second)
      if studyName(f) != s.name || date1 != date2 
        valid = false
        @info "file $path(e) is not valid"
        @show studyName(f), s.name, date1, date2 
      #else
      #  @info "file $path(e) is valid"
      end
    end
  end
  return valid
end

function enforceStudy(s::Study)
  changeParam(s, "/study/name", s.name) 
  changeParam(s, "/study/time", string(s.date)) 
  return
end
