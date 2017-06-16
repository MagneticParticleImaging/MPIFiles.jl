include("Jcampdx.jl")

export BrukerFile, BrukerFileMeas, BrukerFileCalib, latin1toutf8

function latin1toutf8(str::AbstractString)
  buff = Char[]
  for c in str.data
    push!(buff,c)
  end
  string(buff...)
end

function latin1toutf8(str::Void)
  println(stacktrace())
end

abstract BrukerFile <: MPIFile

type BrukerFileMeas <: BrukerFile
  path::String
  params::JcampdxFile
  paramsProc::JcampdxFile
  methodRead
  acqpRead
  visupars_globalRead
  recoRead
  methrecoRead
  visuparsRead
  maxEntriesAcqp
end

type BrukerFileCalib <: BrukerFile
  path::String
  params::JcampdxFile
  paramsProc::JcampdxFile
  methodRead
  acqpRead
  visupars_globalRead
  recoRead
  methrecoRead
  visuparsRead
  maxEntriesAcqp
end

_iscalib(path::String) = isfile(joinpath(path,"pdata", "1", "systemMatrix"))

function (::Type{BrukerFile})(path::String; isCalib=_iscalib(path), maxEntriesAcqp=2000)
  params = JcampdxFile()
  paramsProc = JcampdxFile()

  if isCalib
    return BrukerFileCalib(path, params, paramsProc, false, false, false,
               false, false, false, maxEntriesAcqp)
  else
    return BrukerFileMeas(path, params, paramsProc, false, false, false,
               false, false, false, maxEntriesAcqp)
  end
end

function (::Type{BrukerFile})()
  params = JcampdxFile()
  paramsProc = JcampdxFile()
  return BrukerFileMeas("", params, paramsProc, false, false, false,
             false, false, false, 1)
end

BrukerFileFast(path) = BrukerFile(path, maxEntriesAcqp=400)

function getindex(b::BrukerFile, parameter)
  if !b.acqpRead && ( parameter=="NA" || parameter[1:3] == "ACQ" )
    acqppath = joinpath(b.path, "acqp")
    read(b.params, acqppath, maxEntries=b.maxEntriesAcqp)
    b.acqpRead = true
  elseif !b.methodRead && length(parameter) >= 3 &&
         (parameter[1:3] == "PVM" || parameter[1:3] == "MPI")
    methodpath = joinpath(b.path, "method")
    read(b.params, methodpath)
    b.methodRead = true
  elseif !b.visupars_globalRead && length(parameter) >= 4 &&
         parameter[1:4] == "Visu"
    visupath = joinpath(b.path, "visu_pars")
    read(b.params, visupath, maxEntries=55)
    b.visupars_globalRead = true
  end

  if haskey(b.params, parameter)
    return b.params[parameter]
  else
    return nothing
  end
end

function getindex(b::BrukerFile, parameter, procno::Int64)
  if !b.recoRead && lowercase( parameter[1:4] ) == "reco"
    recopath = joinpath(b.path, "pdata", string(procno), "reco")
    read(b.paramsProc, acqppath, maxEntries=13)
    b.recoRead = true
  elseif !b.methrecoRead && parameter[1:3] == "PVM"
    methrecopath = joinpath(b.path, "pdata", string(procno), "methreco")
    read(b.paramsProc, methrecopath)
    b.methrecoRead = true
  elseif !b.visuparsRead && parameter[1:4] == "Visu"
    visuparspath = joinpath(b.path, "pdata", string(procno), "visu_pars")
    read(b.paramsProc, visuparspath)
    b.visuparsRead = true
  end

  return b.paramsProc[parameter]
end

function Base.show(io::IO, b::BrukerFile)
  print(io, "BrukerFile: ", b.path)
end

# Helper
activeChannels(b::BrukerFile) = [parse(Int64,s) for s=b["PVM_MPI_ActiveChannels"]]
selectedChannels(b::BrukerFile) = b["PVM_MPI_ChannelSelect"] .== "Yes"
selectedReceivers(b::BrukerFile) = b["ACQ_ReceiverSelect"] .== "Yes"

# general parameters
version(b::BrukerFile) = nothing
uuid(b::BrukerFile) = nothing #str2uuid(b["VisuUid"])
time(b::BrukerFile) = nothing

# study parameters
studyName(b::BrukerFile) = string(experimentSubject(b),"_",
                                  latin1toutf8(b["VisuStudyId"]),"_",
                                  b["VisuStudyNumber"])
studyNumber(b::BrukerFile) = parse(Int64,b["VisuStudyNumber"])
studyUuid(b::BrukerFile) = nothing #str2uuid(b["VisuUid"])
studyDescription(b::BrukerFile) = "n.a."

# study parameters
experimentName(b::BrukerFile) = latin1toutf8(b["ACQ_scan_name"])
experimentNumber(b::BrukerFile) = parse(Int64,b["VisuExperimentNumber"])
experimentUuid(b::BrukerFile) = nothing #str2uuid(b["VisuUid"])
experimentDescription(b::BrukerFile) = latin1toutf8(b["ACQ_scan_name"])
experimentSubject(b::BrukerFile) = latin1toutf8(b["VisuSubjectName"])
experimentIsSimulation(b::BrukerFile) = false
experimentIsCalibration(b::BrukerFile) = b["PVM_Matrix"] != nothing
experimentHasProcessing(b::BrukerFile) = experimentIsCalibration(b)
experimentHasReconstruction(b::BrukerFile) = false # fixme later
experimentHasMeasurement(b::BrukerFile) = true

# tracer parameters
tracerName(b::BrukerFile) = [b["PVM_MPI_Tracer"]]
tracerBatch(b::BrukerFile) = [b["PVM_MPI_TracerBatch"]]
tracerVolume(b::BrukerFile) = [parse(Float64,b["PVM_MPI_TracerVolume"])*1e-6]
tracerConcentration(b::BrukerFile) = [parse(Float64,b["PVM_MPI_TracerConcentration"])]
tracerSolute(b::BrukerFile) = ["Fe"]
function tracerInjectionTime(b::BrukerFile)
  initialFrames = b["MPI_InitialFrames"]
  if initialFrames == nothing
    return [acqStartTime(b)]
  else
    return [acqStartTime(b) + Dates.Millisecond(
       round(Int64,parse(Int64, initialFrames)*dfPeriod(b)*1000 ) )]
  end
end
tracerVendor(b::BrukerFile) = ["n.a."]

# scanner parameters
scannerFacility(b::BrukerFile) = latin1toutf8(b["ACQ_institution"])
scannerOperator(b::BrukerFile) = latin1toutf8(b["ACQ_operator"])
scannerManufacturer(b::BrukerFile) = "Bruker/Philips"
scannerModel(b::BrukerFile) = b["ACQ_station"]
scannerTopology(b::BrukerFile) = "FFP"

# acquisition parameters
function acqStartTime(b::BrukerFile)
  acq = b["ACQ_time"] #b["VisuAcqDate"]
  DateTime( replace(acq[2:search(acq,'+')-1],",",".") )
end
acqNumFrames(b::BrukerFile) = Int64(b["ACQ_jobs"][1][8])
acqFramePeriod(b::BrukerFile) = dfPeriod(b) * acqNumAverages(b)
acqNumPatches(b::BrukerFile) = 1
acqNumPeriods(b::BrukerFile) = 1
acqNumAverages(b::BrukerFile) = parse(Int,b["NA"])

function acqNumBGFrames(b::BrukerFile)
  n = b["PVM_MPI_NrBackgroundMeasurementCalibrationAllScans"]
  if n == nothing
    return 0
  else
    return parse(Int64,n)
  end
end
acqGradient(b::BrukerFile) = addTrailingSingleton([-0.5, -0.5, 1.0].*
      parse(Float64,b["ACQ_MPI_selection_field_gradient"]),2)
function acqOffsetField(b::BrukerFile) #TODO NOT correct
  voltage = [parse(Float64,s) for s in b["ACQ_MPI_frame_list"]]
  calibFac = [2.5/49.45, -2.5*0.008/-22.73, 2.5*0.008/-22.73, 1.5*0.0094/13.2963]
  return addTrailingSingleton( Float64[voltage[d]*calibFac[d] for d=2:4],2)
end
acqOffsetFieldShift(b::BrukerFile) = acqOffsetField(b) ./ acqGradient(b)


# drive-field parameters
dfNumChannels(b::BrukerFile) = sum( selectedReceivers(b)[1:3] .== true )
   #sum( dfStrength(b)[1,:,1] .> 0) #TODO Not sure about this
dfStrength(b::BrukerFile) = addTrailingSingleton( addLeadingSingleton(
  [parse(Float64,s) for s = b["ACQ_MPI_drive_field_strength"] ] *1e-3, 2), 3)
dfPhase(b::BrukerFile) = dfStrength(b) .*0 .+  1.5707963267948966 # Bruker specific!
dfBaseFrequency(b::BrukerFile) = 2.5e6
dfCustomWaveform(b::BrukerFile) = nothing
dfDivider(b::BrukerFile) = addTrailingSingleton([102; 96; 99],2)
dfWaveform(b::BrukerFile) = "sine"
dfPeriod(b::BrukerFile) = parse(Float64,b["PVM_MPI_DriveFieldCycle"]) / 1000
# The following takes faked 1D/2D measurements into account
#function dfPeriod(b::BrukerFile)
#  df = dfStrength(b)
#  return lcm(  dfDivider(b)[ (df .>= 0.0000001) .* selectedChannels(b) ] ) / 2.5e6  # in ms!
#end


# receiver parameters
rxNumChannels(b::BrukerFile) = sum( selectedReceivers(b)[1:3] .== true )
rxBandwidth(b::BrukerFile) = parse(Float64,b["PVM_MPI_Bandwidth"])*1e6
rxNumSamplingPoints(b::BrukerFile) = parse(Int64,b["ACQ_size"][1])
rxTransferFunction(b::BrukerFile) = nothing

# measurements
measUnit(b::BrukerFile) = "a.u."
measDataConversionFactor(b::BrukerFileMeas) = [1.0/acqNumAverages(b), 0.0]
measDataConversionFactor(b::BrukerFileCalib) = [1.0, 0.0]

function measData(b::BrukerFileMeas, frames=1:acqNumFrames(b), patches=1:acqNumPatches(b),
                  receivers=1:rxNumChannels(b))

  dataFilename = joinpath(b.path,"rawdata")
  dType = acqNumAverages(b) == 1 ? Int16 : Int32

  raw = Rawfile(dataFilename, dType,
             [rxNumSamplingPoints(b),rxNumChannels(b),acqNumFrames(b)],
             extRaw=".job0") #Int or Uint?
  data = raw[:,receivers,frames]

  return reshape(data,size(data,1),size(data,2),1,size(data,3))
end

function measData(b::BrukerFileCalib, frames=1:acqNumFrames(b), patches=1:acqNumPatches(b),
                  receivers=1:rxNumChannels(b))

  sfFilename = joinpath(b.path,"pdata", "1", "systemMatrix")
  nFreq = rxNumFrequencies(b)

  data = Rawfile(sfFilename, Complex128,
               [prod(calibSize(b)),nFreq,rxNumChannels(b)], extRaw="")
  S = data[]
  scale!(S,1.0/acqNumAverages(b))
  S = reshape(S,size(S,1),size(S,2),size(S,3),1)

  bgFilename = joinpath(b.path,"pdata", "1", "background")

  bgdata = Rawfile(bgFilename, Complex128,
               [acqNumBGFrames(b),nFreq,rxNumChannels(b)], extRaw="")[]
  scale!(bgdata,1.0/acqNumAverages(b))
  #bgdata = permutedims(bgdata,[3,1,2])
  bgdata = reshape(bgdata,size(bgdata,1),size(bgdata,2),size(bgdata,3),1)
  return cat(1,S,bgdata)
end

systemMatrixWithBG(b::BrukerFileCalib) = measData(b)

function systemMatrix(b::BrukerFileCalib, rows, bgCorrection=true)

  localSFFilename = bgCorrection ? "systemMatrixBG" : "systemMatrix"
  sfFilename = joinpath(b.path,"pdata", "1", localSFFilename)
  nFreq = rxNumFrequencies(b)

  data = Rawfile(sfFilename, Complex128,
                 [prod(calibSize(b)),nFreq*rxNumChannels(b)], extRaw="")
  S = data[:,rows]
  scale!(S,1.0/acqNumAverages(b))
  return S
end

measIsFourierTransformed(b::BrukerFileMeas) = false
measIsFourierTransformed(b::BrukerFileCalib) = true
measIsTFCorrected(b::BrukerFile) = false
measIsBGCorrected(b::BrukerFileMeas) = false
# We have it, but by default we pretend that it is not applied
measIsBGCorrected(b::BrukerFileCalib) = false

measIsTransposed(b::BrukerFileMeas) = false
measIsTransposed(b::BrukerFileCalib) = true

measIsFramePermutation(b::BrukerFileMeas) = false
measIsFramePermutation(b::BrukerFileCalib) = true

function measIsBGFrame(b::BrukerFileMeas)
  if !experimentIsCalibration(b)
    return zeros(Bool, acqNumFrames(b))
  else
    isBG = zeros(Bool, acqNumFrames(b))
    increment = parse(Int,b["PVM_MPI_BackgroundMeasurementCalibrationIncrement"])+1
    isBG[1:increment:end] = true
    return isBG
  end
end

# We assume here that the BG frames are at the end
measIsBGFrame(b::BrukerFileCalib) =
   cat(1,zeros(Bool,acqNumFGFrames(b)),ones(Bool,acqNumBGFrames(b)))

measFramePermutation(b::BrukerFileMeas) = nothing
function measFramePermutation(b::BrukerFileCalib)
  perm1=cat(1,measFGFrameIdx(b),measBGFrameIdx(b))
  perm2=cat(1,fgFramePermutation(b),(length(perm1)-acqNumBGFrames(b)+1):length(perm1))
  permJoint = perm1[perm2]
  return permJoint
end

function fgFramePermutation(b::BrukerFile)
  N = calibSize(b)
  counter = 1
  idx = zeros(Int,N...)
  for z=1:N[3]
    for y=1:N[2]
      y_ = mod(z,2)==0 ? N[2]-y+1 : y
      for x=1:N[1]
        x_ = mod(y,2)==0 ? N[1]-x+1 : x
        idx[x_,y_,z] = counter
        counter += 1
      end
    end
  end
  return vec(idx)
end


measIsSpectralLeakageCorrected(b::BrukerFile) = get(b.params, "ACQ_MPI_spectral_cleaningl", "No") != "No"
measIsFrequencySelection(b::BrukerFile) = false

# calibrations
function calibSNR(b::BrukerFile)
  snrFilename = joinpath(b.path,"pdata", "1", "snr")
  data = Rawfile(snrFilename, Float64, [rxNumFrequencies(b),rxNumChannels(b)], extRaw="")
  return addTrailingSingleton(data[],3)
end
calibFov(b::BrukerFile) = [parse(Float64,s) for s = b["PVM_Fov"] ] * 1e-3
calibFovCenter(b::BrukerFile) =
          [parse(Float64,s) for s = b["PVM_MPI_FovCenter"] ] * 1e-3
calibSize(b::BrukerFile) = [parse(Int64,s) for s in b["PVM_Matrix"]]
calibOrder(b::BrukerFile) = "xyz"
calibPositions(b::BrukerFile) = nothing
calibOffsetField(b::BrukerFile) = nothing
calibDeltaSampleSize(b::BrukerFile) = nothing #TODO
calibMethod(b::BrukerFile) = "robot"


# additional functions that should be implemented by an MPIFile
filepath(b::BrukerFile) = b.path



# special additional methods


function sfPath(b::BrukerFile)
  tmp = b["PVM_MPI_FilenameSystemMatrix",1]
  tmp[1:search(tmp,"/pdata")[1]]
end

### The following is for field measurements from Alex Webers method
numCurrentSettings(b::BrukerFile) = parse(Int64,b["MPI_NrCurrentSettings"])
function currentSetting(b::BrukerFile)
  c = Float64[]
  for s in b["MPI_CurrentSetting"]
    append!(c,s)
  end
  return reshape(c,4,div(length(c),4))
end
ballRadius(b::BrukerFile) = parse(Float64,b["MPI_BallRadius"])
numLatitude(b::BrukerFile) = parse(Int64,b["MPI_NrLatitude"])
numMeridian(b::BrukerFile) = parse(Int64,b["MPI_NrMeridian"])
