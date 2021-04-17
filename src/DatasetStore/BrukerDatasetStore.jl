struct BrukerDatasetStore <: DatasetStore
  path::String
end

const BrukerStore = BrukerDatasetStore("/opt/mpidata")

path(e::Experiment{BrukerDatasetStore}) = joinpath( path(e.study), string(e.num) )
path(s::Study{BrukerDatasetStore}, numExp::Integer) = joinpath(path(s),string(numExp))
iscalib(e::Experiment{BrukerDatasetStore}) = _iscalib(path(e))

###  Implementations of abstract interfaces ###

readonly(::BrukerDatasetStore) = true

studydir(d::BrukerDatasetStore) = d.path
calibdir(d::BrukerDatasetStore) = error("BrukerDatasetStore has no calibdir")

function getStudy(d::BrukerDatasetStore, studyfolder::String)
  study = nothing
  studypath = joinpath(d.path,studyfolder)
  if !ishidden(studypath) && isdir(studypath)
    w = split(studyfolder,'_')
    if length(w) >= 5 && length(w[1])==8 # only these can be study folders
      # w_ = w[1:end-2]
      # date = w[1]
      # date = string(date[1:4],"/",date[5:6],"/",date[7:8])

      w = split(studyfolder,'_')
      dateStr = w[1]
      timeStr = w[2]
      date = DateTime(string(dateStr[1:4],"-",dateStr[5:6],"-",dateStr[7:8],"T",
			   timeStr[1:2],":",timeStr[3:4],":",timeStr[5:6]))

      j = JcampdxFile()
      subjfile = string(studypath,"/subject")
      if isfile(subjfile)
        read(j,string(studypath,"/subject"),maxEntries=14) #magic number...
        name = latin1toutf8(j["SUBJECT_study_name"])
        # name = string(latin1toutf8(j["SUBJECT_name_string"]),
        #              "_",latin1toutf8(j["SUBJECT_study_name"]),
        #              "_",latin1toutf8(j["SUBJECT_study_nr"]))
        s1 = latin1toutf8(j["SUBJECT_id"])
	      s2 = latin1toutf8(j["SUBJECT_name_string"])
        subject = (s1 == s2) ? s1 : s1*s2
      else
        # Workaround if no subject file is present => use first dataset
        # and derive the study from the Brukerfile
        r = readdir(studypath)

        found = false
        for file in r

          if tryparse(Int64,file) != nothing
            b = BrukerFileFast(joinpath(studypath, file ))
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
      study = Study(d, name, studyfolder, date, subject)
    end
  end
  return study
end

function getExperiments(s::Study{BrukerDatasetStore})
  files = findBrukerFiles(path(s)) # make me fast

  experiments = Experiment[]

  for file in files
    prefix, ext = splitext(splitpath(file)[end])
    num = tryparse(Int64, prefix)
    if num != nothing 
      #try
      exp = getExperiment(s,num)
      push!(experiments, exp)
      #catch e
      #  @debug "" e
      #end
    end
  end
  return experiments
end

@static if Sys.isunix()
  function findBrukerFiles(path::AbstractString, mindepth::Int=1, maxdepth::Int=2)
    candidatePaths = split(read(`find $path -maxdepth $maxdepth -mindepth $mindepth -type d`,String),"\n")[1:end-1]
    mask = zeros(Bool,length(candidatePaths))
    for (i,candidatePath) in enumerate(candidatePaths)
      if isfile(joinpath(candidatePath,"acqp")) &&
         isfile(joinpath(candidatePath,"method")) &&
         isfile(joinpath(candidatePath,"visu_pars"))
        mask[i] = true
      end
    end
    return String.(candidatePaths[mask])
  end
else
  function findBrukerFiles(path::AbstractString)
    files = readdir(path)
    bfiles = String[]
    for file in files
      if isdir(joinpath(path,file))
       try
        if isfile(joinpath(path,file,"acqp")) &&
           isfile(joinpath(candidatePath,"method")) &&
           isfile(joinpath(candidatePath,"visu_pars"))
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
  return bfiles
  end
end

