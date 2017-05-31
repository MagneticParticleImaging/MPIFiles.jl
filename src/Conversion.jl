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
  params["studyNumber"] = studyNumber(f)
  params["studyDescription"] = studyDescription(f)

  # experiment parameters
  params["experimentName"] = experimentName(f)
  params["experimentNumber"] = experimentNumber(f)
  params["experimentDescription"] = experimentDescription(f)
  params["experimentSubject"] = experimentSubject(f)
  params["experimentIsSimulation"] = experimentIsSimulation(f)
  params["experimentIsCalibration"] = experimentIsCalibration(f)

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
  setparam!(params, "acqOffsetFieldShift", acqOffsetFieldShift(f))

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
  params["measDataConversionFactor"] = measDataConversionFactor(f)
  setparam!(params, "measData", measData(f))
  setparam!(params, "measIsBG", measIsBG(f))

  if params["experimentIsCalibration"]
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
  write(file, "/study/number", get(params,"studyNumber",0))
  write(file, "/study/description", get(params,"studyDescription","n.a."))

  # experiment parameters
  write(file, "/experiment/name", get(params,"experimentName","default") )
  write(file, "/experiment/number", get(params,"experimentNumber",0))
  write(file, "/experiment/description", get(params,"experimentDescription","n.a."))
  write(file, "/experiment/subject", get(params,"experimentSubject","n.a."))
  write(file, "/experiment/isSimulation", Int(get(params,"experimentIsSimulation",false)))

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
  if haskey(params,"acqOffsetFieldShift")
    write(file, "/acquisition/offsetFieldShift", params["acqOffsetFieldShift"])
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
  write(file, "/measurement/dataConversionFactor",  params["measDataConversionFactor"])
  if haskey(params,"measData")
    write(file, "/measurement/data", params["measData"])
  end
  if haskey(params,"measIsBG")
    write(file, "/measurement/isBackgroundData",  convert(Array{Int8},params["measIsBG"]))
  end

  # calibrations
  if params["experimentIsCalibration"]
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
