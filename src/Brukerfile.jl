include("Jcampdx.jl")

export BrukerFile, BrukerFileMeas, BrukerFileCalib, BrukerFileFast, latin1toutf8

function latin1toutf8(str::AbstractString)
  buff = Char[]
  for c in Vector{UInt8}(str)
    push!(buff,c)
  end
  string(buff...)
end

function latin1toutf8(str::Void)
  println(stacktrace())
end

@compat abstract type BrukerFile <: MPIFile end

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

function getindex(b::BrukerFile, parameter)#::String
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
    if isfile(visupath)
      read(b.params, visupath, maxEntries=55)
      b.visupars_globalRead = true
    end
  end

  if haskey(b.params, parameter)
    return b.params[parameter]
  else
    return ""
  end
end

function getindex(b::BrukerFile, parameter, procno::Int64)::String
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
    if isfile(visuparspath)
      read(b.paramsProc, visuparspath)
      b.visuparsRead = true
    end
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
experimentIsCalibration(b::BrukerFile) = _iscalib(b.path)
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
  if initialFrames == ""
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
scannerName(b::BrukerFile) = b["ACQ_station"]
scannerTopology(b::BrukerFile) = "FFP"

# acquisition parameters
function acqStartTime(b::BrukerFile)
  acq = b["ACQ_time"] #b["VisuAcqDate"]
  DateTime( replace(acq[2:search(acq,'+')-1],",",".") )
end
function acqNumFrames(b::BrukerFileMeas)
  M = Int64(b["ACQ_jobs"][1][8])
  return div(M,acqNumPeriods(b))
end
function acqNumFrames(b::BrukerFileCalib)
  M = parse(Int64,b["PVM_MPI_NrCalibrationScans"])
  A = parse(Int64,b["PVM_MPI_NrBackgroundMeasurementCalibrationAdditionalScans"])
  return div(M-A,acqNumPeriods(b))
end
acqFramePeriod(b::BrukerFile) = dfPeriod(b) * acqNumAverages(b)
function _acqNumPatches(b::BrukerFile)
  M = b["MPI_NSteps"]
  return (M == "") ? 1 : parse(Int64,M)
end
function acqNumPeriods(b::BrukerFile)
  M = b["MPI_RepetitionsPerStep"]
  N = _acqNumPatches(b)
  return (M == "") ? N : N*parse(Int64,M)
end

acqNumAverages(b::BrukerFile) = parse(Int,b["NA"])

function acqNumBGFrames(b::BrukerFile)
  n = b["PVM_MPI_NrBackgroundMeasurementCalibrationAllScans"]
  a = b["PVM_MPI_NrBackgroundMeasurementCalibrationAdditionalScans"]
  if n == ""
    return 0
  else
    return parse(Int64,n)-parse(Int64,a)
  end
end
acqGradient(b::BrukerFile) = addTrailingSingleton([-0.5, -0.5, 1.0].*
      parse(Float64,b["ACQ_MPI_selection_field_gradient"]),2)

function acqOffsetField(b::BrukerFile) #TODO NOT correct
  if b["MPI_FocusFieldX"] == ""
    voltage = [parse(Float64,s) for s in b["ACQ_MPI_frame_list"]]
    voltage = reshape(voltage,4,:)
    voltage = repeat(voltage,inner=(1,div(acqNumPeriods(b),_acqNumPatches(b))))
    calibFac = [2.5/49.45, 0.5*(-2.5)*0.008/-22.73, 0.5*2.5*0.008/-22.73, 1.5*0.0094/13.2963]
    return Float64[voltage[d,j]*calibFac[d] for d=2:4, j=1:acqNumPeriods(b)]
  else
    return repeat(1e-3*cat(2,[-parse(Float64,a) for a in b["MPI_FocusFieldX"]],
                 [-parse(Float64,a) for a in b["MPI_FocusFieldY"]],
                 [-parse(Float64,a) for a in b["MPI_FocusFieldZ"]])',inner=(1,div(acqNumPeriods(b),_acqNumPatches(b))))
  end
end


# drive-field parameters
dfNumChannels(b::BrukerFile) = sum( selectedReceivers(b)[1:3] .== true )
   #sum( dfStrength(b)[1,:,1] .> 0) #TODO Not sure about this
dfStrength(b::BrukerFile) = repeat(addTrailingSingleton( addLeadingSingleton(
  [parse(Float64,s) for s = b["ACQ_MPI_drive_field_strength"] ] *1e-3, 2), 3),
                    inner=(1,1,acqNumPeriods(b)))
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
rxInductionFactor(b::BrukerFile) = nothing
rxUnit(b::BrukerFile) = "a.u."
rxDataConversionFactor(b::BrukerFileMeas) =
                 repeat([1.0/acqNumAverages(b), 0.0], outer=(1,rxNumChannels(b)))
rxDataConversionFactor(b::BrukerFileCalib) =
                 repeat([1.0, 0.0], outer=(1,rxNumChannels(b)))

function measData(b::BrukerFileMeas, frames=1:acqNumFrames(b), periods=1:acqNumPeriods(b),
                  receivers=1:rxNumChannels(b))

  dataFilename = joinpath(b.path,"rawdata.job0")
  dType = acqNumAverages(b) == 1 ? Int16 : Int32

  s = open(dataFilename)
  raw = Mmap.mmap(s, Array{dType,4},
             (rxNumSamplingPoints(b),rxNumChannels(b),acqNumPeriods(b),acqNumFrames(b)))
  data = raw[:,receivers,periods,frames]
  close(s)

  return reshape(data, rxNumSamplingPoints(b), length(receivers),length(periods),length(frames))
end

function measData(b::BrukerFileCalib, frames=1:acqNumFrames(b), periods=1:acqNumPeriods(b),
                  receivers=1:rxNumChannels(b))

  sfFilename = joinpath(b.path,"pdata", "1", "systemMatrix")
  nFreq = rxNumFrequencies(b)

  s = open(sfFilename)
  data = Mmap.mmap(s, Array{Complex128,4}, (prod(calibSize(b)),nFreq,rxNumChannels(b),1))
  #S = data[:,:,:,:]
  S = map(Complex64, data)
  close(s)
  scale!(S,1.0/acqNumAverages(b))

  bgFilename = joinpath(b.path,"pdata", "1", "background")

  s = open(bgFilename)
  data = Mmap.mmap(s, Array{Complex128,4}, (acqNumBGFrames(b),nFreq,rxNumChannels(b),1))
  #bgdata = data[:,:,:,:]
  bgdata = map(Complex64, data)
  close(s)
  scale!(bgdata,1.0/acqNumAverages(b))
  return cat(1,S,bgdata)
end

systemMatrixWithBG(b::BrukerFileCalib) = measData(b)

function systemMatrix(b::BrukerFileCalib, rows, bgCorrection=true)

  localSFFilename = bgCorrection ? "systemMatrixBG" : "systemMatrix"
  sfFilename = joinpath(b.path,"pdata", "1", localSFFilename)
  nFreq = rxNumFrequencies(b)

  s = open(sfFilename)
  data = Mmap.mmap(s, Array{Complex128,2}, (prod(calibSize(b)),nFreq*rxNumChannels(b)))
  S = data[:,rows]
  close(s)
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
    # If the file is not a calibration file we cannot say if any particular scans
    # were BG scans
    return zeros(Bool, acqNumFrames(b))
  else
    # In case of a calibration file we know the particular indices corresponding
    # to BG measurements
    isBG = zeros(Bool, acqNumFrames(b))
    increment = parse(Int,b["PVM_MPI_BackgroundMeasurementCalibrationIncrement"])+1
    isBG[1:increment:end] = true

    return isBG
  end
end

# If the file is considered to be a calibration file, we will load
# the measurement in a processed form. In that case the BG measurements
# will be put at the end of the frame dimension.
measIsBGFrame(b::BrukerFileCalib) =
   cat(1,zeros(Bool,acqNumFGFrames(b)),ones(Bool,acqNumBGFrames(b)))

# measurements are not permuted
measFramePermutation(b::BrukerFileMeas) = nothing
# calibration scans are permuted
function measFramePermutation(b::BrukerFileCalib)
  # The following is a trick to obtain the permutation applied to the measurements
  # in a calibration measurement.
  bMeas = BrukerFile(b.path, isCalib=false)

  perm1=cat(1,measFGFrameIdx(bMeas),measBGFrameIdx(bMeas))
  perm2=cat(1,fgFramePermutation(bMeas),(length(perm1)-acqNumBGFrames(bMeas)+1):length(perm1))
  permJoint = perm1[perm2]
  return permJoint
end

# TODO the following requires a test
function fgFramePermutation(b::BrukerFile)
  N = tuple(calibSize(b)...)

  perm = Array{Int}(N)
  for i in CartesianRange(N)
    idx = [i[k] for k=1:length(i)]
    for d=2:3
      if isodd(sum(idx[d:3])-length(idx[d:3]))
        idx[d-1] = N[d-1] + 1 - idx[d-1]
      end
    end
    perm[i] = sub2ind(N,idx...)
  end
  return vec(perm)
end


measIsSpectralLeakageCorrected(b::BrukerFile) = get(b.params, "ACQ_MPI_spectral_cleaningl", "No") != "No"
measIsFrequencySelection(b::BrukerFile) = false

# calibrations
function calibSNR(b::BrukerFile)
  snrFilename = joinpath(b.path,"pdata", "1", "snr")
  s = open(snrFilename)
  data = Mmap.mmap(s, Array{Float64,3}, (rxNumFrequencies(b),rxNumChannels(b),1))
  snr = data[:,:,:]
  close(s)

  return snr
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
