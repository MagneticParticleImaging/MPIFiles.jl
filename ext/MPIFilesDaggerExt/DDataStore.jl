struct DaggerDatasetStore{T, C <: Dagger.Chunk{T}} <: DDatasetStore{T}
  chunk::Dagger.Chunk{T}
  worker::Int64
  # Differentiate between bruker and mdf store?
  function MPIFiles.DDatasetStore(args...; worker::Int64)
    chunk = Dagger.@mutable worker = worker MDFDatasetStore(args...)
    return new{MDFDatasetStore, typeof(chunk)}(chunk, worker)
  end
end
MPIFiles.worker(store::DaggerDatasetStore) = store.worker


function changeStore(study::Study, store::DatasetStore)
  return Study(store, study.name, study.foldername, study.date, study.subject, study.uuid)
end
function changeStore(exp::Experiment, store::DatasetStore)
  return Experiment(changeStore(exp.study, store), exp.num, exp.name, exp.numFrames, exp.df, exp.sfGradient, exp.numAverages, exp.operator, exp.time)
end

function MPIFiles.path(e::Experiment{<:DaggerDatasetStore}) 
  return fetch(Dagger.spawn(getMDFStore(e).chunk) do store
    exp = changeStore(e, store)
    return path(exp)
  end)
end
function MPIFiles.path(s::Study{<:DaggerDatasetStore}, numExp::Integer) 
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    return path(study, numExp)
  end)
end
MPIFiles.readonly(store::DaggerDatasetStore) = fetch(Dagger.spawn(readonly, store.chunk))
MPIFiles.studydir(store::DaggerDatasetStore) = fetch(Dagger.spawn(studydir, store.chunk))
Base.empty!(store::DaggerDatasetStore) = fetch(Dagger.spawn(empty!, store.chunk))


function MPIFiles.iscalib(e::Experiment{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(e).chunk) do store
    exp = changeStore(e, store)
    return MPIFiles.iscalib(exp)
  end)
end
function MPIFiles.changeParam(e::Experiment{<:DaggerDatasetStore}, paramName::AbstractString, paramValue)
  wait(Dagger.spawn(getMDFStore(e).chunk) do store
    exp = changeStore(e, store)
    changeParam(path(exp), paramName, paramValue)
  end)
end

MPIFiles.calibdir(store::DaggerDatasetStore) = fetch(Dagger.spawn(calibdir, store.chunk))

function MPIFiles.getMDFStore(study::Study{<:DaggerDatasetStore})
    return study.store
end

function MPIFiles.getMDFStore(experiment::Experiment{<:DaggerDatasetStore})
    return getMDFStore(experiment.study)
end

function MPIFiles.getStudy(d::DaggerDatasetStore, studyfolder::String)
  return fetch(Dagger.spawn(d.chunk) do store
    study = getStudy(store, studyfolder)
    return changeStore(study, d)
  end)
end

function MPIFiles.getCalibStudy(d::DaggerDatasetStore)
  return fetch(Dagger.spawn(d.chunk) do store
    study = getCalibStudy(store)
    return changeStore(study, d)
  end)
end

function MPIFiles.getExperiments(s::Study{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    exps = getExperiments(study)
    return map(exp -> changeStore(exp, getMDFStore(s)), exps)
  end)
end

function MPIFiles.getMDFStudyFolderName(s::Study{<:DaggerDatasetStore}) 
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    return getMDFStudyFolderName(study)
  end)
end

function MPIFiles.addStudy(d::DaggerDatasetStore, s::Study{<:DaggerDatasetStore})
  # Fetch to get errors
  return fetch(Dagger.spawn(d.chunk) do store
    study = changeStore(s, store)
    addStudy(store, study)
    return nothing
  end)
end

MPIFiles.getNewCalibNum(d::DaggerDatasetStore) = fetch(Dagger.spawn(getNewCalibNum, d.chunk))
MPIFiles.getNewCalibPath(d::DaggerDatasetStore) = fetch(Dagger.spawn(MPIFiles.getNewCalibPath, d.chunk))