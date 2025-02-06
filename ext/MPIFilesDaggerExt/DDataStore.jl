struct DaggerDatasetStore{T, C <: Dagger.Chunk{T}} <: DDatasetStore{T}
  chunk::Dagger.Chunk{T}
  worker::Int64
end
# Differentiate between bruker and mdf store?
function DDataStore(args...; worker::Int64)
  chunk = Dagger.@mutable worker = worker MDFDatasetStore(args...)
  return DaggerDataStore(chunk, worker)
end
MPIFiles.worker(store::DaggerDataStore) = store.worker


function changeStore(study::Study{DaggerDatasetStore}, store::DatasetStore)
  return Study(store, study.name, study.foldername, study.date, study.subject, study.uuid)
end
function changeStore(exp::Experiment{DaggerDatasetStore}, store::DatasetStore)
  return Experiment(changeStore(exp.study, store), exp.num, exp.name, exp.numFrames, exp.df, exp.sfGradient, exp.numAverages, exp.operator, exp.time)
end

function path(e::Experiment{DaggerDataStore}) 
  return fetch(Dagger.spawn(getMDFStore(e)) do store
    exp = changeStore(e, store)
    return path(exp)
  end)
end
function path(s::Study{DaggerDataStore}, numExp::Integer) 
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    return path(study, numExp)
  end)
end
MPIFiles.readonly(store::DaggerDataStore) = fetch(Dagger.@spawn readonly(store.chunk))
MPIFiles.studydir(store::DaggerDataStore) = fetch(Dagger.@spawn studydir(store.chunk))
Base.empty!(store::DaggerDataStore) = fetch(Dagger.@spawn empty!(store.chunk))


function iscalib(e::Experiment{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(e)) do store
    exp = changeStore(e, store)
    return iscalib(e)
  end)
end
# TODO
changeParam(e::Experiment{DaggerDataStore}, paramName::AbstractString, paramValue) =
               changeParam(path(e), paramName, paramValue)

studydir(d::DaggerDataStore) = fetch(Dagger.@spawn studydir(store.chunk))
calibdir(d::DaggerDataStore) = fetch(Dagger.@spawn calibdir(store.chunk))

function getMDFStore(study::Study{DaggerDataStore})
    return study.store
end

function getMDFStore(experiment::Experiment{DaggerDataStore})
    return getMDFStore(experiment.study)
end

function getStudy(d::DaggerDataStore, studyfolder::String)
  return fetch(Dagger.spawn(d.chunk) do store
    study = getStudy(store, studyfolder)
    return changeStore(study, d)
  end)
end

function getCalibStudy(d::DaggerDataStore)
  return fetch(Dagger.spawn(d.chunk) do store
    study = getCalibStudy(store)
    return changeStore(study, d)
  end)
end

function getExperiments(s::Study{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    exps = getExperiments(s)
    return map(exp -> changeStore(exp, s), exps)
  end)
end

getMDFStudyFolderName(study::Study{DaggerDatasetStore}) = fetch(Dagger.@spawn getMDFStudyFolderName(study.chunk))

function addStudy(d::DaggerDataStore, s::Study{DaggerDatasetStore})
  # Fetch to get errors
  return fetch(Dagger.spawn(d.chunk) do store
    study = changeStore(s, store)
    addStudy(store, study)
    return nothing
  end)
end

MPIFiles.getNewCalibNum(d::DaggerDataStore) = fetch(Dagger.@spawn getNewCalibNum(d.chunk))
MPIFiles.getNewCalibPath(d::DaggerDataStore) = fetch(Dagger.@spawn getNewCalibPath(d.chunk))