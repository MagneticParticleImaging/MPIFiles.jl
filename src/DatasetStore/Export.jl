export exportData, makeTarGzip, createArtifact

function exportData(e::Experiment, mdf::MDFDatasetStore, study::Union{Nothing,Study}=nothing; logging=false, 
                    storeForwardRef=false, keepExpNum=false, kargs...)
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
    if study == nothing
      name = studyName(f)
      subject = experimentSubject(f)
      date = studyTime(f)
      s = Study(mdf, name; date=date, subject=subject)
    else
      s = study
    end
    if keepExpNum  
      expNum = experimentNumber(f) # or e.num
      exportpath = joinpath(path(s),string(expNum)*".mdf")
    else
      expNum = getNewExperimentNum(s)
      exportpath = getNewExperimentPath(s)
    end
    saveasMDF(exportpath, f; experimentNumber = expNum, kargs...)
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
      write(io,"$(path(e)), $exportpath\n")
    end
  end

  return exportpath
end

function exportData(s::Study, mdf::MDFDatasetStore, args...; 
                    nums::Vector{Int}=zeros(Int,0),
                    kargs...)
  exps = getExperiments(s)
  for e in exps
    # nums indicates which experiments are exported
    (length(nums) == 0 || e.num in nums) && exportData(e, mdf, args...; kargs...)
  end

  return nothing
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

function createArtifact(d::DatasetStore, url="https://")
  makeTarGzip(d)
  tarball = joinpath(splitpath(d.path)...)*".tar.gz" 
  artifact_toml = joinpath(d.path, "..", "Artifacts.toml")
  tarball_hash = bytes2hex(GitTools.blob_hash(tarball))
  name = splitpath(d.path)[end]

  hash = create_artifact() do artifact_dir
    unpack(tarball, artifact_dir)
  end
  bind_artifact!(artifact_toml, name, hash;
                 download_info=[(url, tarball_hash)],lazy=true, force=true)
end

function validate(d::DatasetStore)
  return all([validate(s) for s in getStudies(d)])
end
