export Study, getStudy, getStudies, validate, enforceStudy, changeStudy

function path(s::Study{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    return path(study)
  end)
end

function getStudies(d::DaggerDataStore, name::String="")
  return fetch(Dagger.spawn(d.chunk) do store
    studies = getStudies(d, name)
    return map(s -> changeStore(s, d), studies)
  end)
end

function remove(s::Study{DaggerDataStore})
  fetcH(Dagger.spawn(getMDFStore(s).chunk) do store
    study = changeStore(s, store)
    remove(study)
    return nothing
  end)
end


function getNewExperimentNum(s::Study{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    getNewExperimentNum(store)
  end)
end

function getNewExperimentPath(s::Study{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    getNewExperimentPath(store)
  end)
end

function changeStudy(s::Study{DaggerDataStore}, newStudyName::AbstractString; date::DateTime = now()) 
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    changeStudy(store, newStudyName; date)
  end)
end

function validate(s::Study{DaggerDataStore})
  return fetch(Dagger.spawn(getMDFStore(s)) do store
    study = changeStore(s, store)
    validate(store)
  end)
end

