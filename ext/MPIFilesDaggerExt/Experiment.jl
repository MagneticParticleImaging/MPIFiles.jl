function MPIFile(e::Experiment{DaggerDatasetStore}; kwargs...)
  # Worker last to overwrite potential worker in kwargs
  return DMPIFile(path(e); kwargs..., worker = worker(e.study.store))
end

function getExperiment(study::Study{DaggerDataStore}, numExp::Integer)
  return fetch(Dagger.spawn(getMDFStore(study)) do store
    s = changeStore(study, store)
    exp = getExperiment(s, numExp)
    return changeStore(exp, getMDFStore(study))
  end)
end

function remove(exp::Experiment)
  if isfile(path(exp))
    rm(path(exp))
  end
end