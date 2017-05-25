module MPIFiles

using Graphics: @mustimplement

import Base: ndims, call, time

### export list ###

export MPIFile

# general parameters
export version, uuid

# study parameters
export studyName, studyExperiment, studyDescription, studySubject,
       studyIsSimulation, studyIsCalibration

# tracer parameters
export tracerName, tracerBatch, tracerVolume, tracerConcentration,
       tracerSolute, tracerTime

# scanner parameters
export scannerFacility, scannerOperator, scannerManufacturer, scannerModel,
       scannerTopology

# acquisition parameters
export acqTime, acqNumFrames, acqNumBGFrames, acqFramePeriod, acqNumPatches,
       acqGradient, acqOffsetField, acqFov, acqFovCenter

# drive-field parameters
export dfNumChannels, dfStrength, dfPhase, dfBaseFrequency, dfCustomWaveform,
       dfDivider, dfWaveform, dfPeriod

# receiver parameters
export rxNumChannels, rxNumAverages, rxBandwidth, rxNumSamplingPoints, rxFrequencies,
       rxTransferFunction

# measurements
export measUnit, measRawDataConversion,
       measData, measDataTimeOrder, measBGData, measBGDataTimeOrder

# calibrations
export calibSystemMatrixData, calibSNR, calibFov, calibFovCenter, calibSize,
       calibOrder, calibPositions, calibOffsetField, calibDeltaSampleSize,
       calibMethod

# reconstruction results
export recoData, recoFov, recoFovCenter, recoSize, recoOrder, recoPositions

### Interface of an MPIFile ###

abstract MPIFile

# general parameters
@mustimplement version(f::MPIFile)
@mustimplement uuid(f::MPIFile)
@mustimplement time(f::MPIFile)

# study parameters
@mustimplement studyName(f::MPIFile)
@mustimplement studyExperiment(f::MPIFile)
@mustimplement studyDescription(f::MPIFile)
@mustimplement studySubject(f::MPIFile)
@mustimplement studyIsSimulation(f::MPIFile)
@mustimplement studyIsCalibration(f::MPIFile)

# tracer parameters
@mustimplement tracerName(f::MPIFile)
@mustimplement tracerBatch(f::MPIFile)
@mustimplement tracerVolume(f::MPIFile)
@mustimplement tracerConcentration(f::MPIFile)
@mustimplement tracerSolute(f::MPIFile)
@mustimplement tracerTime(f::MPIFile)

# scanner parameters
@mustimplement scannerFacility(f::MPIFile)
@mustimplement scannerOperator(f::MPIFile)
@mustimplement scannerManufacturer(f::MPIFile)
@mustimplement scannerModel(f::MPIFile)
@mustimplement scannerTopology(f::MPIFile)

# acquisition parameters
@mustimplement acqTime(f::MPIFile)
@mustimplement acqNumFrames(f::MPIFile)
@mustimplement acqNumBGFrames(f::MPIFile)
@mustimplement acqFramePeriod(f::MPIFile)
@mustimplement acqNumPatches(f::MPIFile)
@mustimplement acqGradient(f::MPIFile)
@mustimplement acqOffsetField(f::MPIFile)
@mustimplement acqFov(f::MPIFile)
@mustimplement acqFovCenter(f::MPIFile)

# drive-field parameters
@mustimplement dfNumChannels(f::MPIFile)
@mustimplement dfStrength(f::MPIFile)
@mustimplement dfPhase(f::MPIFile)
@mustimplement dfBaseFrequency(f::MPIFile)
@mustimplement dfCustomWaveform(f::MPIFile)
@mustimplement dfDivider(f::MPIFile)
@mustimplement dfWaveform(f::MPIFile)
@mustimplement dfPeriod(f::MPIFile)

# receiver properties
@mustimplement rxNumChannels(f::MPIFile)
@mustimplement rxNumAverages(f::MPIFile)
@mustimplement rxBandwidth(f::MPIFile)
@mustimplement rxNumSamplingPoints(f::MPIFile)
@mustimplement rxFrequencies(f::MPIFile)
@mustimplement rxTransferFunction(f::MPIFile)

# measurements
@mustimplement measUnit(f::MPIFile)
@mustimplement measRawDataConversion(f::MPIFile)
@mustimplement measData(f::MPIFile)
@mustimplement measDataTimeOrder(f::MPIFile)
@mustimplement measBGData(f::MPIFile)
@mustimplement measBGDataTimeOrder(f::MPIFile)

# calibrations
@mustimplement calibSystemMatrixData(f::MPIFile)
@mustimplement calibSNR(f::MPIFile)
@mustimplement calibFov(f::MPIFile)
@mustimplement calibFovCenter(f::MPIFile)
@mustimplement calibSize(f::MPIFile)
@mustimplement calibOrder(f::MPIFile)
@mustimplement calibPositions(f::MPIFile)
@mustimplement calibOffsetField(f::MPIFile)
@mustimplement calibDeltaSampleSize(f::MPIFile)
@mustimplement calibMethod(f::MPIFile)




# reconstruction results
@mustimplement recoData(f::MPIFile)
@mustimplement recoFov(f::MPIFile)
@mustimplement recoFovCenter(f::MPIFile)
@mustimplement recoSize(f::MPIFile)
@mustimplement recoOrder(f::MPIFile)
@mustimplement recoPositions(f::MPIFile)














#TODO
@mustimplement filepath(f::MPIFile)
@mustimplement getTimeDataHandle(f::MPIFile)

# non abstract functions
export ndims, voxelSize, voxelVolume, consistenceCheck, loadHeader

### These are generic and can also be used for HDF5 MPI file
ndims(b::MPIFile) = sum( (dfStrength(b) .> 0.00001) )
voxelSize(b::MPIFile) = fov(b) ./ gridSize(b)
voxelVolume(b::MPIFile) = prod( voxelSize(b) ) * 1000 #in Liter
dfFov(b::MPIFile) = 2*dfStrength(b) ./ abs( sfGradient(b) ) # MPI_DFFov???
measPath{T<:MPIFile}(bs::Vector{T}) = [measPath(b) for b in bs]
gridSize{T<:MPIFile}(bs::Vector{T}) = [gridSize(b) for b in bs]
gridSizeCommon{T<:MPIFile}(bs::Vector{T}) = gridSize(bs[1])
gridSizeCommon(bs::MPIFile) = gridSize(bs)
numReceivers{T<:MPIFile}(bs::Vector{T}) = minimum([numReceivers(b) for b in bs])
sfGradient(b::MPIFile, d::Integer) = sfGradient(b::MPIFile)[d]
hasSpectralCleaning(b::MPIFile) = false


function consistenceCheck(bSF::MPIFile, bMeas::MPIFile)
  gSF = sfGradient(bSF,3)
  gMeas = sfGradient(bMeas,3)
  if gSF != gMeas
    warn("The gradient strength of the system matrix ($gSF T/m) does not fit to the measurements ($gMeas T/m)!")
  end

  dfSF = dfStrength(bSF)
  dfMeas = dfStrength(bMeas)
  if dfSF != dfMeas
    warn("The drive-field strength of the system matrix ($dfSF mT) does not fit to the measurements ($dfMeas mT)!")
  end

end

function consistenceCheck{T<:MPIFile}(bSFs::Vector{T}, bMeas::MPIFile)
  for bSF in bSFs
    consistenceCheck(bSF,bMeas)
  end
end

function consistenceCheck{T<:MPIFile}(bSF::MPIFile, bMeass::Vector{T})
  for bMeas in bMeass
    consistenceCheck(bSF,bMeas)
  end
end
function consistenceCheck{T<:MPIFile}(bSFs::Vector{T}, bMeass::Vector{T})
  for i = 1:length(bMeass)
    bMeas=bMeass[i]
    bSF=bSFs[i]
    consistenceCheck(bSF,bMeas)
  end
end

# Don't want to make this public and hide this therefore
stepsize(f::MPIFile) = 1

function loadHeader(b::MPIFile)
  header = Dict{String,Any}()
  try
    # study properties
    header["subjectName"] = subjectName(b)
    header["studyName"] = studyName(b)
    header["measPath"] = measPath(b)
    header["expno"] = expno(b)
    header["description"] = description(b)
    header["isReference"] = isReference(b)
    header["isSimulation"] = isSimulation(b)

    # scanner properties
    header["facility"] = facility(b)
    header["operator"] = operator(b)
    header["manufacturer"] = manufacturer(b)
    header["model"] = model(b)
    header["topology"] = topology(b)

    # tracer properties
    header["tracer"] = tracer(b)
    header["tracerVendor"] = tracerVendor(b)
    header["tracerBatch"] = tracerBatch(b)
    header["tracerVolume"] = tracerVolume(b)
    header["tracerConcentration"] = tracerConcentration(b)
    header["timeTracerInjection"] = timeTracerInjection(b)

    # general acquisition properties
    header["sfGradient"] = sfGradient(b)
    header["numScans"] = numScans(b)
    header["acqDate"] = acqDate(b)
    header["numPatches"] = numPatches(b)
    header["framePeriod"] = framePeriod(b)

    # drive-field properties
    header["dfStrength"] = dfStrength(b)
    header["dfFov"] = dfFov(b)
    header["dfcycle"] = dfcycle(b)
    header["dfBaseFrequency"] = dfBaseFrequency(b)
    header["dfDivider"] = dfDivider(b)
    header["dfRepetitionTime"] = dfRepetitionTime(b)
    header["ffPos"] = ffPos(b)
    header["numAverages"] = numAverages(b)
    header["numDfChannels"] = numDfChannels(b)


    # receiver properties
    header["bandwidth"] = bandwidth(b)
    #header["frequencies"] = frequencies(b) #This is not really needed later
    header["numFreq"] = numFreq(b)
    header["numReceivers"] = numReceivers(b)
    header["numTimePoints"] = numTimePoints(b)


    # TODO: move the following to Analyze???
    dateStr, timeStr = split("$(acqDate(b))","T")
    dateStr = prod(split(dateStr,"-"))
    timeStr = split(timeStr,".")[1]
    timeStr = prod(split(timeStr,":"))

    header["date"] = dateStr
    header["time"] = timeStr
  catch err
   println(err)
  end
  return header
end

# Generate image header that can be used to promote
# general properties of a measurement to a reconstruction
function generateHeaderDict(bSF::MPIFile, b::MPIFile)

  #header["spatialorder"] TODO

  header = loadHeader(b)

  header["datatype"] = "MPI"
  header["dim"] = ndims(b)
  header["size"] = gridSize(bSF)
  #header["pixelspacing"] = voxelSize(bSF) ./ sfGradient(b,3) .* sfGradient(bSF,3)
  #header["offset"] = ffPos(b)

  header
end

# multi patch fallback
generateHeaderDict(bSF::MPIFile, b::Vector) = generateHeaderDict(bSF,b[1])

### Concrete implementations ###
include("Brukerfile.jl")
include("MDF.jl")
include("Image.jl")
include("Conversion.jl")
include("ConversionIMT.jl")

# This dispatches on the file extension and automatically
# generates the correct type
@compat function (::Type{MPIFile})(filename::AbstractString)
  filenamebase, ext = splitext(filename)
  if ext == ".mdf" || ext == ".hdf" || ext == ".h5"
    return MDFFile(filename)
  else
    return BrukerFile(filename)
  end
end

@compat function (::Type{MPIFile})(filenames::Vector)
  return map(x->MPIFile(x),filenames)
end

include("Background.jl")
include("Measurements.jl")
include("SystemMatrix.jl")
include("FrequencyFilter.jl")
include("Reconstruction.jl")
include("TrustedFOV.jl")


### Misc functions ###
include("Misc.jl")


end # module
