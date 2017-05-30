# This file contains routines to generate MDF files
export saveasMDF, loadFullDataset

function setparam!(params::Dict, parameter, value)
  if value != nothing
    params[parameter] = value
  end
end

function loadFullDataset(f)
  params = Dict{String,Any}()

  # general parameters
  params["version"] = version(f)
  setparam!(params, "uuid", uuid(f))
  setparam!(params, "time", time(f))

  # study parameters
  params["studyName"] = studyName(f)
  params["studyExperiment"] = studyExperiment(f)
  params["studyDescription"] = studyDescription(f)
  params["studySubject"] = studySubject(f)
  params["studyIsSimulation"] = studyIsSimulation(f)
  params["studyIsCalibration"] = studyIsCalibration(f)

  # tracer parameters
  params["tracerName"] = tracerName(f)
  params["tracerBatch"] = tracerBatch(f)
  params["tracerVendor"] = tracerVendor(f)
  params["tracerVolume"] = tracerVolume(f)
  params["tracerConcentration"] = tracerConcentration(f)
  params["tracerSolute"] = tracerSolute(f)
  params["tracerInjectionTime"] = tracerInjectionTime(f)

  # scanner parameters
  params["scannerFacility"] = scannerFacility(f)
  params["scannerOperator"] = scannerOperator(f)
  params["scannerManufacturer"] = scannerManufacturer(f)
  params["scannerModel"] = scannerModel(f)
  params["scannerTopology"] = scannerTopology(f)

  # acquisition parameters
  params["acqNumFrames"] = acqNumFrames(f)
  params["acqNumBGFrames"] = acqNumBGFrames(f)
  params["acqFramePeriod"] = acqFramePeriod(f)
  params["acqNumPatches"] = acqNumPatches(f)
  params["acqStartTime"] = acqStartTime(f)
  setparam!(params, "acqGradient", acqGradient(f))
  setparam!(params, "acqOffsetField", acqOffsetField(f))
  setparam!(params, "acqFov", acqFov(f))
  setparam!(params, "acqFovCenter", acqFovCenter(f))

  # drivefield parameters
  params["dfNumChannels"] = dfNumChannels(f)
  params["dfStrength"] = dfStrength(f)
  params["dfPhase"] = dfPhase(f)
  params["dfBaseFrequency"] = dfBaseFrequency(f)
  params["dfDivider"] = dfDivider(f)
  params["dfPeriod"] = dfPeriod(f)
  params["dfWaveform"] = dfWaveform(f)
  if params["dfWaveform"] == "custom"
    params["dfCustomWaveform"] = dfCustomWaveform(f)
  end

  # receiver parameters
  params["rxNumChannels"] = rxNumChannels(f)
  params["rxNumAverages"] = rxNumAverages(f)
  params["rxBandwidth"] = rxBandwidth(f)
  params["rxNumSamplingPoints"] = rxNumSamplingPoints(f)
  setparam!(params, "rxTransferFunction", rxTransferFunction(f))

  # measurement
  params["measUnit"] = measUnit(f)
  params["measRawDataConversion"] = measRawDataConversion(f)
  setparam!(params, "measData", measData(f))
  setparam!(params, "measDataTimeOrder", measDataTimeOrder(f))
  setparam!(params, "measBGData", measBGData(f))
  setparam!(params, "measBGDataTimeOrder", measBGDataTimeOrder(f))

  if params["studyIsCalibration"]
    setparam!(params, "calibSystemMatrixData", calibSystemMatrixData(f))
    setparam!(params, "calibSNR", calibSNR(f))

    params["calibFov"] = calibFov(f)
    params["calibFovCenter"] = calibFovCenter(f)
    params["calibSize"] = calibSize(f)
    params["calibOrder"] = calibOrder(f)
    setparam!(params, "calibPositions", calibPositions(f))
    setparam!(params, "calibOffsetField", calibOffsetField(f))
    setparam!(params, "calibDeltaSampleSize", calibDeltaSampleSize(f))
    params["calibMethod"] = calibMethod(f)
  end

  return params
end

function saveasMDF(filenameOut::String, filenameIn::String)
  saveasMDF(filenameOut, loadFullDataset(MPIFile(filenameIn)) )
end

function saveasMDF(filename::String, params::Dict)
  h5open(filename, "w") do file
    saveasMDF(file, params)
  end
end

function saveasMDF(file::HDF5File, params::Dict)
  # general parameters
  write(file, "/version", "2.0")
  write(file, "/uuid", get(params,"uuid",hex(rand(UInt128))) )
  write(file, "/time", "$( get(params,"time", Dates.unix2datetime(time())) )")

  # study parameters
  write(file, "/study/name", get(params,"studyName","default") )
  write(file, "/study/experiment", get(params,"studyExperiment",0))
  write(file, "/study/description", get(params,"studyDescription","n.a."))
  write(file, "/study/subject", get(params,"studySubject","n.a."))
  write(file, "/study/isSimulation", Int(get(params,"studyIsSimulation",false)))
  write(file, "/study/isCalibration", Int(haskey(params,"calibSize")))

  # tracer parameters
  write(file, "/tracer/name", get(params,"tracerName","n.a") )
  write(file, "/tracer/batch", get(params,"tracerBatch","n.a") )
  write(file, "/tracer/vendor", get(params,"tracerVendor","n.a") )
  write(file, "/tracer/volume", get(params,"tracerVolume",0.0))
  write(file, "/tracer/concentration", get(params,"tracerConcentration",0.0) )
  write(file, "/tracer/solute", get(params,"tracerSolute","Fe") )
  write(file, "/tracer/injectionTime", "$( get(params,"tracerInjectionTime", Dates.unix2datetime(time())) )")

  # scanner parameters
  write(file, "/scanner/facility", get(params,"scannerFacility","n.a") )
  write(file, "/scanner/operator", get(params,"scannerOperator","n.a") )
  write(file, "/scanner/manufacturer", get(params,"scannerManufacturer","n.a"))
  write(file, "/scanner/model", get(params,"scannerModel","n.a"))
  write(file, "/scanner/topology", get(params,"scannerTopology","FFP"))

  # acquisition parameters
  write(file, "/acquisition/numFrames", get(params,"acqNumFrames",1))
  write(file, "/acquisition/numBackgroundFrames", get(params,"acqNumBGFrames",0))
  write(file, "/acquisition/framePeriod", get(params,"acqFramePeriod",0.0))
  write(file, "/acquisition/numPatches", get(params,"acqNumPatches",1))
  write(file, "/acquisition/startTime", "$( get(params,"acqStartTime", Dates.unix2datetime(time())) )")

  if haskey(params,"acqGradient")
    write(file, "/acquisition/gradient", params["acqGradient"])
  end
  if haskey(params,"acqOffsetField")
    write(file, "/acquisition/offsetField", params["acqOffsetField"])
  end
  if haskey(params,"acqFov")
    write(file, "/acquisition/fieldOfView", params["acqFov"])
  end
  if haskey(params,"acqFovCenter")
    write(file, "/acquisition/fieldOfViewCenter", params["acqFovCenter"])
  end

  # drivefield parameters
  write(file, "/acquisition/drivefield/numChannels", size(params["dfStrength"],2) )
  write(file, "/acquisition/drivefield/strength", params["dfStrength"])
  write(file, "/acquisition/drivefield/phase", params["dfPhase"])
  write(file, "/acquisition/drivefield/baseFrequency", params["dfBaseFrequency"])
  write(file, "/acquisition/drivefield/divider", params["dfDivider"])
  write(file, "/acquisition/drivefield/period", params["dfPeriod"])
  write(file, "/acquisition/drivefield/waveform", params["dfWaveform"])
  if params["dfWaveform"] == "custom"
    write(file, "/acquisition/drivefield/customWaveform", params["dfCustomWaveform"])
  end

  # receiver parameters
  write(file, "/acquisition/receiver/numChannels", params["rxNumChannels"])
  write(file, "/acquisition/receiver/numAverages",  params["rxNumAverages"])
  write(file, "/acquisition/receiver/bandwidth", params["rxBandwidth"])
  write(file, "/acquisition/receiver/numSamplingPoints", params["rxNumSamplingPoints"])

  if haskey(params,"rxTransferFunction")
    tf = params["rxTransferFunction"]
    tfR = reinterpret(Float64, tf, (2, size(tf)...))
    write(file, "/acquisition/receiver/transferFunction", tfR)
  end

  # measurements
  write(file, "/measurement/unit",  params["measUnit"])
  write(file, "/measurement/rawDataConversion",  params["measRawDataConversion"])
  if haskey(params,"measData")
    write(file, "/measurement/data",  params["measData"])
  end
  if haskey(params,"measDataTimeOrder")
    write(file, "/measurement/dataTimeOrder",  params["measDataTimeOrder"])
  end
  if haskey(params,"measBGData")
    write(file, "/measurement/backgroundData",  params["measBGData"])
  end
  if haskey(params,"measBGDataTimeOrder")
    write(file, "/measurement/backgroundDataTimeOrder",  params["measBGDataTimeOrder"])
  end

  # calibrations
  if params["studyIsCalibration"]
    if haskey(params,"calibSystemMatrixData")
      S = params["calibSystemMatrixData"]
      S = reinterpret(typeof((S[1]).re),S,(2,size(S)...))
      write(file, "/calibration/systemMatrixData", S)
    end
    if haskey(params,"calibSNR")
      write(file, "/calibration/snr",  params["calibSNR"])
    end
    write(file, "/calibration/fieldOfView",  params["calibFov"])
    write(file, "/calibration/fieldOfViewCenter",  params["calibFovCenter"])
    write(file, "/calibration/size",  params["calibSize"])
    write(file, "/calibration/order",  params["calibOrder"])
    if haskey(params,"calibPositions")
      write(file, "/calibration/positions",  params["calibPositions"])
    end
    if haskey(params,"calibOffsetField")
      write(file, "/calibration/offsetField",  params["calibOffsetField"])
    end
    if haskey(params,"calibDeltaSampleSize")
      write(file, "/calibration/deltaSampleSize",  params["calibDeltaSampleSize"])
    end
    write(file, "/calibration/method",  params["calibMethod"])
  end


  #TODO reconstruction results

end
