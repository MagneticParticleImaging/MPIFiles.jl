# This file contains routines to generate MDF files
export saveasMDF, loadFullDataset, loadMetadata

function setparam!(params::Dict, parameter, value)
  if value != nothing
    params[parameter] = value
  end
end

function loadFullDataset(f)
  params = loadMetadata(f)

  # call API function and store result in a parameter Dict
  if experimentHasMeasurement(f)
    for op in [:measUnit, :measData, :measDataConversionFactor, :measIsBG]
      setparam!(params, string(op), eval(op)(f))
    end
  end

  if experimentHasProcessing(f)
    for op in [:procData, :procIsFourierTransformed, :procIsTFCorrected,
               :procIsAveraged, :procIsFramesSelected, :procIsBGCorrected,
               :procIsTransposed, :procFramePermutation]
      setparam!(params, string(op), eval(op)(f))
    end
  end

  if experimentIsCalibration(f)
    for op in [:calibSNR, :calibFov, :calibFovCenter,
               :calibSize, :calibOrder, :calibPositions, :calibOffsetField,
               :calibDeltaSampleSize, :calibMethod]
      setparam!(params, string(op), eval(op)(f))
    end
  end

  if experimentHasReconstruction(f)
    for op in [:recoData, :recoSize, :recoFov, :recoFovCenter, :recoOrder,
               :recoPositions, :recoParameters]
      setparam!(params, string(op), eval(op)(f))
    end
  end

  return params
end

function loadMetadata(f)
  params = Dict{String,Any}()

  # call API function and store result in a parameter Dict
  for op in [:version, :uuid, :time, :studyName, :studyNumber, :studyUuid, :studyDescription,
            :experimentName, :experimentNumber, :experimentUuid, :experimentDescription, 
            :experimentSubject,
            :experimentIsSimulation, :experimentIsCalibration, :experimentHasProcessing,
            :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
            :tracerSolute, :tracerInjectionTime,
            :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerModel,
            :scannerTopology, :acqNumFrames, :acqNumBGFrames, :acqFramePeriod,
            :acqNumPatches, :acqStartTime, :acqGradient, :acqOffsetField, :acqOffsetFieldShift,
            :dfNumChannels, :dfStrength, :dfPhase, :dfBaseFrequency, :dfDivider,
            :dfPeriod, :dfWaveform, :rxNumChannels, :rxNumAverages, :rxBandwidth,
            :rxNumSamplingPoints, :rxTransferFunction]
    setparam!(params, string(op), eval(op)(f))
  end


  if params["dfWaveform"] == "custom"
    params["dfCustomWaveform"] = dfCustomWaveform(f)
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

hasKeyAndValue(paramDict,param) = haskey(paramDict, param) && paramDict[param] != nothing

function writeIfAvailable(file, paramOut, paramDict, paramIn )
  if hasKeyAndValue(paramDict, paramIn)
    write(file, paramOut, paramDict[paramIn])
  end
end

function saveasMDF(file::HDF5File, params::Dict)
  # general parameters
  write(file, "/version", "2.0")
  write(file, "/uuid", string(get(params,"uuid",Base.Random.uuid4() )))
  write(file, "/time", "$( get(params,"time", Dates.unix2datetime(time())) )")

  # study parameters
  write(file, "/study/name", get(params,"studyName","default") )
  write(file, "/study/number", get(params,"studyNumber",0))
  if hasKeyAndValue(params,"studyUuid")
    studyUuid = params["studyUuid"]
  else
    studyUuid = Base.Random.uuid4()
  end
  write(file, "/study/uuid", string(studyUuid))
  write(file, "/study/description", get(params,"studyDescription","n.a."))

  # experiment parameters
  write(file, "/experiment/name", get(params,"experimentName","default") )
  write(file, "/experiment/number", get(params,"experimentNumber",0))
  if hasKeyAndValue(params,"experimentUuid")
    expUuid = params["experimentUuid"]
  else
    expUuid = Base.Random.uuid4()
  end
  write(file, "/experiment/uuid", string(expUuid))
  write(file, "/experiment/description", get(params,"experimentDescription","n.a."))
  write(file, "/experiment/subject", get(params,"experimentSubject","n.a."))
  write(file, "/experiment/isSimulation", Int8(get(params,"experimentIsSimulation",false)))

  # tracer parameters
  write(file, "/tracer/name", get(params,"tracerName","n.a") )
  write(file, "/tracer/batch", get(params,"tracerBatch","n.a") )
  write(file, "/tracer/vendor", get(params,"tracerVendor","n.a") )
  write(file, "/tracer/volume", get(params,"tracerVolume",0.0))
  write(file, "/tracer/concentration", get(params,"tracerConcentration",0.0) )
  write(file, "/tracer/solute", get(params,"tracerSolute","Fe") )
  tr = [string(t) for t in get(params,"tracerInjectionTime", [Dates.unix2datetime(time())]) ]
  write(file, "/tracer/injectionTime", tr)

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

  writeIfAvailable(file, "/acquisition/gradient", params, "acqGradient")
  writeIfAvailable(file, "/acquisition/offsetField", params, "acqOffsetField")
  writeIfAvailable(file, "/acquisition/offsetFieldShift", params, "acqOffsetFieldShift")

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

  if hasKeyAndValue(params,"rxTransferFunction")
    tf = params["rxTransferFunction"]
    tfR = reinterpret(Float64, tf, (2, size(tf)...))
    write(file, "/acquisition/receiver/transferFunction", tfR)
  end

  # measurements
  if hasKeyAndValue(params, "measData")
    write(file, "/measurement/unit",  params["measUnit"])
    write(file, "/measurement/dataConversionFactor",  params["measDataConversionFactor"])
    write(file, "/measurement/data", params["measData"])
    if hasKeyAndValue(params,"measIsBG")
      write(file, "/measurement/isBackgroundData",  convert(Array{Int8},params["measIsBG"]))
    end
  end

  # processing
  if params["experimentHasProcessing"]
    if params["procIsTransposed"]
      S = params["procData"]
      S = reinterpret(typeof((S[1]).re),S,(2,size(S)...))
      write(file, "/processing/data", S)
    else
      write(file, "/processing/data", S)
    end
    write(file, "/processing/isFourierTransformed", Int8(params["procIsFourierTransformed"]))
    write(file, "/processing/isTransferFunctionCorrected", Int8(params["procIsTFCorrected"]))
    write(file, "/processing/isAveraged",  Int8(params["procIsAveraged"]))
    write(file, "/processing/isFramesSelected", Int8(params["procIsFramesSelected"]))
    write(file, "/processing/isBackgroundCorrected",  Int8(params["procIsBGCorrected"]))
    write(file, "/processing/isTransposed",  Int8(params["procIsTransposed"]))
    if haskey(params,"procFramePermutation")
      write(file, "/processing/framePermutation",  params["procFramePermutation"])
    end
  end

  # calibrations
  if params["experimentIsCalibration"]
    if hasKeyAndValue(params,"calibSNR")
      write(file, "/calibration/snr",  params["calibSNR"])
    end
    write(file, "/calibration/fieldOfView",  params["calibFov"])
    write(file, "/calibration/fieldOfViewCenter",  params["calibFovCenter"])
    write(file, "/calibration/size",  params["calibSize"])
    write(file, "/calibration/order",  params["calibOrder"])
    if hasKeyAndValue(params,"calibPositions")
      write(file, "/calibration/positions",  params["calibPositions"])
    end
    if hasKeyAndValue(params,"calibOffsetField")
      write(file, "/calibration/offsetField",  params["calibOffsetField"])
    end
    if hasKeyAndValue(params,"calibDeltaSampleSize")
      write(file, "/calibration/deltaSampleSize",  params["calibDeltaSampleSize"])
    end
    write(file, "/calibration/method",  params["calibMethod"])
  end

  # reconstruction
  if hasKeyAndValue(params, "recoData")
    write(file, "/reconstruction/data", params["recoData"])
    write(file, "/reconstruction/fieldOfView", params["recoFov"])
    write(file, "/reconstruction/fieldOfViewCenter", params["recoFovCenter"])
    write(file, "/reconstruction/size", params["recoSize"])
    write(file, "/reconstruction/order", get(params,"recoOrder", "xyz"))
    if hasKeyAndValue(params,"recoPositions")
      write(file, "/reconstruction/positions", params["recoPositions"])
    end
    if hasKeyAndValue(params,"recoParameters")
      saveParams(file, "/reconstruction/parameters", params["recoParameters"])
    end
  end

end
