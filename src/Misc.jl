

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
