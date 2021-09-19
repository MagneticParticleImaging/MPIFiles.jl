export Experiment, getExperiment, getExperiments, getNewExperimentNum, getNewExperimentPath, getStudy

# might not be so clever to use an explicit type fields here
# maybe better a dict
struct Experiment{D<:DatasetStore}
  study::Study{D}
  num::Int64
  name::String
  numFrames::Int64
  df::Vector{Float64}
  sfGradient::Float64
  numAverages::Int64
  operator::String
  time::String
  # more ...
end

@mustimplement path(e::Experiment{D}) where D<:DatasetStore

MPIFile(e::Experiment; kargs...) = MPIFile(path(e); kargs...)

function getStudy(experiment::Experiment)
    return experiment.study
end

function getExperiment(s::Study, numExp::Integer)

  path_ = path(s, numExp)

  if !ispath(path_)
    return nothing
  end

  exp = MPIFile(path_, fastMode=true) do b 
    prefix, ext = splitext(path_)
    
    Experiment(s, parse(Int64,last(splitdir(prefix))), #why parse experiment number here????, 
                      string(experimentName(b)), acqNumFrames(b),
                      round.(1000 .* vec(dfStrength(b)[1,:,1]),digits=2), maximum(abs.(acqGradient(b))),
                      acqNumAverages(b), scannerOperator(b), string(acqStartTime(b)))
  end
  return exp
end

function remove(exp::Experiment)
  if isfile(path(exp))
    rm(path(exp))
  end
end