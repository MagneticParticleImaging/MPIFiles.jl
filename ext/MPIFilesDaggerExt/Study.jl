export Study, getStudy, getStudies, validate, enforceStudy, changeStudy

function MPIFiles.path(s::Study{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    return path(study)
  end)
end

function MPIFiles.getStudies(d::DaggerDatasetStore, name::String="")
  return fetch(Dagger.spawn(d.chunk) do store
    studies = getStudies(store, name)
    return map(s -> changeStore(s, d), studies)
  end)
end

function MPIFiles.remove(s::Study{<:DaggerDatasetStore})
  fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    remove(study)
    return nothing
  end)
end


function MPIFiles.getNewExperimentNum(s::Study{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    getNewExperimentNum(study)
  end)
end

function MPIFiles.getNewExperimentPath(s::Study{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    getNewExperimentPath(study)
  end)
end

function MPIFiles.changeStudy(s::Study{<:DaggerDatasetStore}, newStudyName::AbstractString; kwargs...) 
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    changeStudy(study, newStudyName; kwargs...)
  end)
end

function MPIFiles.validate(s::Study{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    validate(store)
  end)
end

