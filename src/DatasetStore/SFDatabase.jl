export generateSFDatabase, loadSFDatabase

function findSFFiles(d::BrukerDatasetStore)
  studies = readdir(d.path)
  bfiles = String[]

  for study in studies
    studypath = joinpath(d.path,study)
    if isdir(studypath)
      experiments = readdir(studypath)
      for exp in experiments
        path = joinpath(d.path,study,exp)
        if _iscalib(path)
          push!(bfiles, path)
        end
      end
    end
  end
  BrukerMDFSFs = readdir("/opt/data/MDF_SFs/")
  for BrukerMDFSF in BrukerMDFSFs
    prefix, ext = splitext(BrukerMDFSF)
     if ext == ".mdf"
      push!(bfiles,joinpath("/opt/data/MDF_SFs/",BrukerMDFSF))
    end
  end
  return bfiles
end

function findSFFiles(d::MDFDatasetStore)
  bfiles = String[]

  path = joinpath(d.path,"calibrations/")

  files = readdir(path)

  for file in files
    prefix, ext = splitext(file)
    if !isdir(file) && tryparse(Int64,prefix) != nothing &&
       (ext == ".mdf" || ext == ".hdf" || ext == ".h5") && !occursin("td.mdf",file)
      try
        push!(bfiles, joinpath(path,file))
      catch e
        @debug "" e
      end
    end
  end

  return bfiles
end



####

function generateSFDatabase(d::DatasetStore, filename::AbstractString)
  fileList = findSFFiles(d)
  A = generateSFDatabase(fileList)
  writedlm(filename, A, ',')
end

function generateSFDatabase(fileList::Vector)

  A = Array{Any}(undef,length(fileList)+1,17)

  # Headerrow
  A[1,1] = "Name"
  A[1,2] = "Gradient"
  A[1,3] = "DFx"
  A[1,4] = "DFy"
  A[1,5] = "DFz"
  A[1,6] = "Size x"
  A[1,7] = "Size y"
  A[1,8] = "Size z"
  A[1,9] = "Bandwidth"
  A[1,10] = "Tracer"
  A[1,11] = "TracerBatch"
  A[1,12] = "DeltaSampleConcentration"
  A[1,13] = "DeltaSampleVolume"
  A[1,14] = "Path"
  A[1,15] = "StartDate"
  A[1,16] = "MeasurementTime"
  A[1,17] = "ExperimentNumber"

  for (k,sf) in enumerate(fileList)
    i=k+1
    b = MPIFile(sf)
    _innerGenerateSFDatabase(A,i,sf,b)
  end
  return A
end

function _innerGenerateSFDatabase(A,i,sf,b)
  A[i,1] = experimentName(b)
  A[i,2] = maximum(acqGradient(b))
  df = vec(dfStrength(b)).*1e3
  A[i,3:5] .= 0.0
  for l=1:min(length(df),3)
    A[i,l+2] = df[l]
  end
  N = calibSize(b)
  A[i,6] = N[1]
  A[i,7] = N[2]
  A[i,8] = N[3]
  A[i,9] = rxBandwidth(b) / 1e6
  A[i,10] = tracerName(b)[1]
  A[i,11] = tracerBatch(b)[1]
  A[i,12] = 0.0#deltaSampleConcentration(b)
  A[i,13] = 0.0#deltaSampleVolume(b)
  A[i,14] = sf #filepath(b)
  A[i,15] = string(acqStartTime(b))
  A[i,16] = 0.0#b["PVM_ScanTimeStr"]
  A[i,17] = experimentNumber(b)
end

function generateSFDatabase(d::MDFDatasetStore)
  oldfile = joinpath(d.path,"SF_DatabaseOld.csv")
  newfile = joinpath(d.path,"SF_Database.csv")
  generateSFDatabase_(d, oldfile, newfile)
end

# HAAACKKK
function generateSFDatabase(d::BrukerDatasetStore)
  oldfile = "/opt/data/SF_DatabaseOld.csv"
  newfile = "/opt/data/SF_Database.csv"
  generateSFDatabase_(d, oldfile, newfile)
end

function generateSFDatabase_(d::DatasetStore, oldfile, newfile)

  if isfile(newfile)
    if isfile(oldfile)
      cp(newfile, oldfile, force=true)
    else
      cp(newfile, oldfile, force=false)
    end
  end

  generateSFDatabase(d, newfile)
end

function loadSFDatabase(d::BrukerDatasetStore)
  if isfile("/opt/data/SF_Database.csv")
    A = readdlm("/opt/data/SF_Database.csv",',')
    if size(A,2) < 16
      A = readdlm("/opt/data/SF_Database.csv",'\t')
    end
    return A
  else
    return nothing
  end
end

function loadSFDatabase(d::MDFDatasetStore)
  files = readdir(calibdir(d))
  @debug "system function database" files
  mdffiles = files[endswith.(files,".mdf") .& .!endswith.(files,"td.mdf")]
  fileList = calibdir(d).*"/".*mdffiles
  A = generateSFDatabase(fileList)
  return A
end

####
