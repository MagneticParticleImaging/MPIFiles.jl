function MPIFiles.MPIFile(e::Experiment{<:DaggerDatasetStore}; kwargs...)
  # Worker last to overwrite potential worker in kwargs
  return DMPIFile(path(e); kwargs..., worker = MPIFiles.worker(getMDFStore(e)))
end

function MPIFiles.getExperiment(study::Study{<:DaggerDatasetStore}, numExp::Integer)
  return fetch(Dagger.spawn(getMDFStore(study).chunk) do store
    s = changeStore(study, store)
    exp = getExperiment(s, numExp)
    return isnothing(exp) ? nothing : changeStore(exp, getMDFStore(study))
  end)
end

function MPIFiles.remove(e::Experiment{<:DaggerDatasetStore})
  return fetch(Dagger.spawn(getMDFStore(e).chunk) do store
    exp = changeStore(e, store)
    remove(exp)
  end)
end