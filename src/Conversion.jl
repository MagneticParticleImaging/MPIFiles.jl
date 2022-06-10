# This file contains routines to generate MDF files
export saveasMDF, loadDataset, loadMetadata, loadMetadataOnline, setparam!, compressCalibMDF
export getFFdataPerPos, prepareAsMDFSingleMeasurement, convertCustomSF, blockAverage


function setparam!(params::Dict, parameter, value)
  if value != nothing
    params[parameter] = value
  end
end

# we do not support all conversion possibilities
function loadDataset(f::MPIFile; frames=1:acqNumFrames(f), applyCalibPostprocessing=false,
                     numPeriodAverages=1, numPeriodGrouping=1, experimentNumber=nothing, kargs...)
  params = loadMetadata(f)

  if experimentNumber != nothing
    params[:experimentNumber] = experimentNumber
  end

  # call API function and store result in a parameter Dict
  if experimentHasMeasurement(f)
    loadMeasParams(f, params, skipMeasData = true)
    if !applyCalibPostprocessing
      if frames!=1:acqNumFrames(f)
        setparam!(params, :measData, measData(f,frames))
        setparam!(params, :acqNumFrames, length(frames))
        setparam!(params, :measIsBGFrame, measIsBGFrame(f)[frames])
      else
        setparam!(params, :measData, measData(f))
      end
    else
        @info "load measurement data"
        data = getMeasurementsFD(f, false, frames=1:acqNumFrames(f), sortFrames=true,
                       spectralLeakageCorrection=false, transposed=true, tfCorrection=false,
                       numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)

          @info size(data)
        setparam!(params, :measData, data)
        setparam!(params, :measIsFourierTransformed, true)
        setparam!(params, :measIsFastFrameAxis, true)
        setparam!(params, :measIsFramePermutation, true)
        setparam!(params, :measFramePermutation, fullFramePermutation(f))

        if numPeriodAverages > 1
          params[:acqNumAverages] *= numPeriodAverages
          params[:acqNumPeriodsPerFrame] = div(params[:acqNumPeriodsPerFrame],numPeriodAverages)
        end
        if numPeriodGrouping > 1
          params[:acqNumPeriodsPerFrame] = div(params[:acqNumPeriodsPerFrame],numPeriodGrouping)
          params[:dfCycle] *= numPeriodGrouping # Not sure about this one

          params[:rxNumSamplingPoints] *= numPeriodGrouping
        end

        setparam!(params, :measIsBGFrame,
          cat(zeros(Bool,acqNumFGFrames(f)),ones(Bool,acqNumBGFrames(f)), dims=1))

        snr = calibSNR(f)
	      if snr == nothing
          @info "calculate SNR"
          snr = calculateSystemMatrixSNR(f, data, numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)
        end
        setparam!(params, :calibSNR, snr)
    end
  end

  loadCalibParams(f, params)
  loadRecoParams(f, params)

  return params
end

const defaultParams =[:version, :uuid, :time, :dfStrength, :acqGradient, :studyName,
          :studyNumber, :studyUuid, :studyTime, :studyDescription,
          :experimentName, :experimentNumber, :experimentUuid, :experimentDescription,
          :experimentSubject,
          :experimentIsSimulation, :experimentIsCalibration,
          :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
          :tracerSolute, :tracerInjectionTime,
          :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName,
          :scannerTopology, :acqNumPeriodsPerFrame, :acqNumAverages,
          :acqStartTime, :acqOffsetField, :acqNumFrames,
          :dfNumChannels, :dfPhase, :dfBaseFrequency, :dfDivider,
          :dfCycle, :dfWaveform, :rxNumChannels, :rxBandwidth,
          :rxNumSamplingPoints, :rxTransferFunction, :rxTransferFunctionFileName, :rxInductionFactor,
          :rxUnit, :rxDataConversionFactor]

function loadMetadata(f, inputParams = MPIFiles.defaultParams)
  params = Dict{Symbol,Any}()
  # call API function and store result in a parameter Dict
  for op in inputParams
    setparam!(params, op, eval(op)(f))
  end
  return params
end

function loadRecoParams(f, params = Dict{Symbol,Any}())
  if experimentHasReconstruction(f)
    for op in [:recoData, :recoSize, :recoFov, :recoFovCenter, :recoOrder,
           :recoPositions, :recoParameters]
      setparam!(params, op, eval(op)(f))
    end
  end

  return params
end

function loadCalibParams(f, params = Dict{Symbol,Any}())
  if experimentIsCalibration(f)
    for op in [:calibFov, :calibFovCenter,
               :calibSize, :calibOrder, :calibPositions, :calibOffsetField,
             :calibDeltaSampleSize, :calibMethod]
      setparam!(params, op, eval(op)(f))
    end
    if !haskey(params, :calibSNR)
      setparam!(params, :calibSNR, calibSNR(f))
    end
  end
  return params
end

function loadMeasParams(f, params = Dict{Symbol,Any}(); skipMeasData = false)
  if experimentHasMeasurement(f)
    for op in [:measIsFourierTransformed, :measIsTFCorrected,
                 :measIsBGCorrected,
                 :measIsFastFrameAxis, :measIsFramePermutation, :measIsFrequencySelection,
                 :measIsSpectralLeakageCorrected,
                 :measFramePermutation, :measIsBGFrame]
      setparam!(params, op, eval(op)(f))
    end
  end

  if !skipMeasData
    setparam!(params, :measData, measData(f))
  end

  if measIsFrequencySelection(f)
    setparam!(params, :measFrequencySelection, measFrequencySelection(f))
  end

  return params
end



function appendBGDataset(params::Dict, filenameBG::String; kargs...)
  params = MPIFile(filenameBG) do fBG 
    appendBGDataset(params, fBG; kargs...)
  end
  return params
end

function appendBGDataset(params::Dict, fBG::MPIFile; frames=1:acqNumFrames(fBG))
  paramsBG = loadDataset(fBG, frames=frames)
  paramsBG[:measIsBGFrame][:] = true

  params[:measData] = cat(4, params[:measData], paramsBG[:measData])
  params[:measIsBGFrame] = cat(params[:measIsBGFrame], paramsBG[:measIsBGFrame], dims=1)
  params[:acqNumFrames] += paramsBG[:acqNumFrames]

  return params
end


isConvertibleToMDF(f::MPIFile) = true
function isConvertibleToMDF(f::BrukerFile)
  # check if raw data is consistent
  if !rawDataLengthConsistent(f::BrukerFile)
    return false
  end
  # check if conversion of method is supported
  #TODO use regex matching if list grows too large
  whitelist = ["User:ukeMPI337", "User:mkaul_MPI", "User:ukeMPINew", "User:ukeMPI335", "User:uke_mpi", "Bruker:MPICalibration", "Bruker:MPI"]
  if !(f["ACQ_method"] in whitelist)
    @warn "conversion of method not supported" f["ACQ_method"]
    return false
  end
  return true
end

function saveasMDF(filenameOut::String, filenameIn::String; kargs...)
  MPIFile(filenameIn) do f
    saveasMDF(filenameOut, f; kargs...)
  end
  return
end

function saveasMDF(filenameOut::String, f::MPIFile; filenameBG = nothing, enforceConversion=false, kargs...)
  
  # This is a hack. Needs to be fixed properly
  if(haskey(kargs, :SNRThresh) || haskey(kargs, :sparsityTrafoRedFactor)) && calibSNR(f) != nothing
    compressCalibMDF(filenameOut, f; kargs...)
    return
  end
  
  if enforceConversion || isConvertibleToMDF(f)
    params = loadDataset(f;kargs...)
    if filenameBG != nothing
      appendBGDataset(params, filenameBG)
    end
    saveasMDF(filenameOut, params)
  else
    error("File not Convertible")
  end
  return
end


function saveasMDF(filename::String, params::Dict)
  # file has to be removed if exists. Otherwise h5create fails.
  isfile(filename) && rm(filename)
  h5open(filename, "w") do file
    saveasMDF(file, params)
  end
end

function compressCalibMDF(filenameOut::String, f::MPIFile; SNRThresh=2.0, kargs...)
  idx = Int64[]

  SNR = calibSNR(f)[:,:,1]
  for k=1:size(SNR,1)
    if maximum(SNR[k,:]) > SNRThresh
      push!(idx, k)
    end
  end

  compressCalibMDF(filenameOut, f, idx; kargs...)
end

function compressCalibMDF(filenamesOut::Vector{String}, f::MultiMPIFile; SNRThresh=2.0, kargs...)
  idx = Int64[]

  # We take the SNR from the first SF
  SNR = calibSNR(f[1])[:,:,1]
  for k=1:size(SNR,1)
    if maximum(SNR[k,:]) > SNRThresh
      push!(idx, k)
    end
  end

  for (i,f_) in enumerate(f)
    compressCalibMDF(filenamesOut[i], f_, idx; kargs...)
  end
end

function compressCalibMDF(filenameOut::String, f::MPIFile, idx::Vector{Int64};
                          sparsityTrafoRedFactor=1.0, sparsityTrafo="DCT-IV", kargs...)
  params = loadMetadata(f)
  loadMeasParams(f, params, skipMeasData = true)
  loadCalibParams(f, params)
  params[:calibSNR] = calibSNR(f)
  loadRecoParams(f, params)

  data = systemMatrixWithBG(f, idx)

  params[:calibSNR] = params[:calibSNR][idx,:,:]
  if haskey(params, :rxTransferFunction)
    params[:rxTransferFunction] = params[:rxTransferFunction][idx,:]
  end
  params[:measIsFrequencySelection] = true
  params[:measFrequencySelection] = idx

  if sparsityTrafoRedFactor == 1.0
    params[:measData] = data
  else
    B = linearOperator(sparsityTrafo, calibSize(f))
    N = prod(calibSize(f))
    NBG = size(data,1) - N
    D = size(data,3)
    P = size(data,4)
    NRed = max(1, floor(Int, sparsityTrafoRedFactor*N))
    dataOut = similar(data, NBG+NRed, length(idx), D, P)
    subsamplingIndices = zeros(Int32, NRed, length(idx), D, P)

    fgdata = data[measFGFrameIdx(f),:,:,:]
    bgdata = data[measBGFrameIdx(f),:,:,:]

    dataOut[(NRed+1):end,:,:,:] = bgdata

    for k=1:length(idx), d=1:D, p=1:P
      I = B * fgdata[:,k,d,p]
      subsamplingIndices[:,k,d,p] = round.(Int32,reverse(sortperm(abs.(I)),dims=1)[1:NRed])
      dataOut[1:NRed,k,d,p] = I[vec(subsamplingIndices[:,k,d,p])]
    end

    params[:measData] = dataOut
    params[:measIsSparsityTransformed] = true
    params[:measSparsityTransformation] = sparsityTrafo
    params[:measSubsamplingIndices] = subsamplingIndices

    bgFrame = zeros(Bool, NRed+NBG)
    bgFrame[(NRed+1):end] .= true
    params[:measIsBGFrame] = bgFrame
    params[:acqNumFrames] = NRed+NBG
  end

  saveasMDF(filenameOut, params)
end


function blockAverage(filenameOut::String, f::MPIFile, numAverages)
	#block averaging is not possible
	if rem(acqNumFrames(f),numAverages) != 0
		throw(DomainError(numAverages, "`numAverages` must be a divider of $(acqNumFrames(f))"))
	end
	
	# is there a way to rule out calibration scans?
	
	# get data, exept for measurement data
	params = loadMetadata(f)
	loadMeasParams(f, params, skipMeasData = true)
	loadCalibParams(f, params)
	loadRecoParams(f, params)

	# uuid and time newly created on saving MDF
	delete!(params,:uuid)
	delete!(params,:time)
	
	# update measurement data
	Nnew = div(acqNumFrames(f),numAverages)
	data = measData(f)
	if measIsSparsityTransformed(f)
		error("No block averaging of sparsity transformed data possible")
	elseif measIsFramePermutation(f)
		error("No block averaging of data with frame permutation possible")
	else
		if measIsFastFrameAxis(f)
			N,J,C,KW = size(data)
			data = reshape(data,numAverages,Nnew,J,C,KW)
			data = mean(data, dims=1) # creates Float64 data if eltype is <: Integer
			data = reshape(data,Nnew,J,C,KW)
		else
			J,C,KW,N = size(data)
			data = reshape(data,J,C,KW,numAverages,Nnew)
			data = mean(data, dims=4) # creates Float64 data if eltype is <: Integer
			data = reshape(data,J,C,KW,Nnew)
		end
		if eltype(data) == Float64
			data = Float32.(data) # Float64 is not necessary
		end
		setparam!(params, :measData, data)
	end
  	
	# update measurement parameters
	p = measIsBGFrame(f)
	p = reshape(p,numAverages,Nnew)
	p = any(p, dims=1) # all blocks with at least one background frames will be labeled as background
	p = reshape(p,Nnew)
	setparam!(params, :measIsBGFrame, p)

	# update acquisition parameters
	setparam!(params, :acqNumFrames, Nnew)
	setparam!(params, :acqNumAverages, acqNumAverages(f)*numAverages)

  return saveasMDF(filenameOut, params)
end

function loadAndProcessFFData(f::BrukerFile, nAverages::Int64, skipSwitchingFrames::Int64;addToEnd=0)
  dataFilename = joinpath(f.path,"rawdata.job0")
  #FileSize = stat(dataFilename).size
  #AllFrames = convert(Int,round(Int,FileSize/26928/4/3))
  #AllFrames = acqNumPeriodsPerPatch(f)*acqNumPeriodsPerFrame(f)
  AllFrames = acqNumPeriodsPerFrame(f)
  (Nx,Ny,Nz) = [length(union(acqOffsetField(f)[ll,1,:])) for ll in collect(1:3)]

  data_ = zeros(ComplexF64,Nx,Ny,Nz,rxNumFrequencies(f),rxNumChannels(f));

#  Pos_=acqOffsetField(f).*[-1,-1,-1]; #Field to position
  Pos=acqOffsetField(f).*[-1,-1,-1]; #Field to position

#  Pos = zeros(3,1,AllFrames);
#  Pos[1,1,:] = kron(vec(Pos_[1,1,:]),ones(2))
#  Pos[2,1,:] = kron(vec(Pos_[2,1,:]),ones(2))
#  Pos[3,1,:] = kron(vec(Pos_[3,1,:]),ones(2))

  dType = Int16 #acqNumAverages(f) == 1 ? Int16 : Int32  # Data are not averaged

  (PosNx,PosNy,PosNz)=[union(Pos[ll,1,:]) for ll in collect(1:3)]

  ds = open(dataFilename)
  if numSubPeriods(f) == 1
     raw = Mmap.mmap(ds, Array{dType,4},(rxNumSamplingPoints(f),1,rxNumChannels(f),AllFrames));
  else
     raw = Mmap.mmap(ds, Array{dType,5},(rxNumSamplingPoints(f),numSubPeriods(f),1,rxNumChannels(f),AllFrames))
     raw = dropdims(sum(raw,dims=2),dims=2)
  end

  for p = collect(1:Nx*Ny*Nz)
#    p = Nx*Ny*(pz-1)+pxy
    st = p*skipSwitchingFrames+(p-1)*(nAverages+addToEnd)+1

    ind = collect(st:st+nAverages-1)

    if any([length(union(Pos[ll,1,ind])) for ll in collect(1:3)].>1)
      error("Averaged data are not at the same FF-Position! Check Parameter Averages and/or skipSwitchingFrames (MPIFiles function loadAndProcessFFData)")
    end
    # sequence is not relevant
    #x = collect(Nx:-1:1)[PosNx.==union(Pos[1,1,ind])]
    #y = collect(Ny:-1:1)[PosNy.==union(Pos[2,1,ind])]
    #z = collect(Nz:-1:1)[PosNz.==union(Pos[3,1,ind])]
    x = collect(1:Nx)[PosNx.==union(Pos[1,1,ind])]
    y = collect(1:Ny)[PosNy.==union(Pos[2,1,ind])]
    z = collect(1:Nz)[PosNz.==union(Pos[3,1,ind])]
    #for ch =1:3
      #dataSC = spectralLeakageCorrectedData(raw[:,1,ch,ind]);
      #data_[x,y,z,:,ch] = rfft(vec(mean(dataSC,dims=2)));
      data_[x,y,z,:,:] = rfft(mean(raw[:,1,:,ind],dims=3),1);
    #end
  end
  data = reshape(data_,:,rxNumFrequencies(f),rxNumChannels(f),1)
  return convert(Array{Complex{Float32},4},data); # allPos x freq x channel x 1
end

function convertCustomSF(filenameOut::String, f::BrukerFile, fBG::BrukerFile,nAverages::Int64,skipSwitchingFrames::Int64; nAveragesBG = nAverages,skipSwitchingFramesBG=skipSwitchingFrames)

  params = loadMetadata(f)
  loadMeasParams(f, params, skipMeasData = true)
  loadCalibParams(f, params)

  paramsBG = loadMetadata(fBG)


  params[:calibSize] = [length(union(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  numFGFrames = prod(params[:calibSize])

  params[:acqGradient] = reshape(params[:acqGradient][:,:,1,1],3,3,1,1)
  params[:dfPhase] =  reshape(params[:dfPhase][:,:,1],1,3,1)
  params[:calibDeltaSampleSize] = [0.0, 0.0, 0.0] #Todo
  #params[:time] = now()
  params[:calibMethod] = "BrukerCustom"
  #params[:version] = v"2.0.0"
  params[:dfStrength] = reshape(params[:dfStrength][:,:,1,1],1,3,1)
  params[:experimentIsCalibration] = true
  #params[:uuid] = uuid4()
  params[:measIsFourierTransformed] = true
  params[:measIsFastFrameAxis] = true
  params[:calibFovCenter] = [mean(extrema(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  params[:acqNumPeriodsPerFrame] = 1

  params[:measIsFramePermutation] = 1

  if rxHasTransferFunction(f) &&!measIsTFCorrected(f)
    params[:transferFunction] = rxTransferFunction(f)
  end


  calibFov=[sum(abs.(extrema(params[:acqOffsetField][ll,1,:])))./abs.(params[:acqGradient][ll,ll]) for ll in collect(1:3)]
  ind=findall(x->x==0.0,calibFov) # if calibFov is 0.0 for 2D Measurements, AxisArray will fail in creating an ImageMeta object
  calibFov[ind].=0.001  # for 1mm
  params[:calibFov] = calibFov
  params[:calibOrder] = "xyz"

  params[:acqOffsetField] = reshape([mean(extrema(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)],3,1,1)
  params[:acqNumAverages] = nAverages

println("Part1")
  # Daten Laden
@time  fgFrames = loadAndProcessFFData(f,nAverages, skipSwitchingFrames) # allePos x freqs x r x 1
@time  bgFramesFull = loadAndProcessFFData(fBG,nAveragesBG, skipSwitchingFramesBG) # allePos x freqs x r x 1
  nAveragesBGHalf = nAveragesBG >1 ? div(nAveragesBG,2) : nAveragesBG
@time  bgFramesHalf = loadAndProcessFFData(fBG,div(nAveragesBG,2), skipSwitchingFramesBG+nAveragesBG-div(nAveragesBG,2)) # allePos x freqs x r x 1

  params[:calibSNR] = calculateSNRCustomSF(f,fgFrames,bgFramesFull,bgFramesHalf)

  (xBG,yBG,zBG) = [length(union(paramsBG[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  (xFG,yFG,zFG) = params[:calibSize]

  bgFullReshaped = reshape(bgFramesFull,xBG,yBG,zBG,rxNumFrequencies(fBG),rxNumChannels(fBG))

  if yBG !=1
     if zBG !=1
       itp = interpolate(bgFullReshaped,(NoInterp(),BSpline(Linear()),BSpline(Linear()),NoInterp(),NoInterp())); #Interpolire bgFramesFull auf
     else
       itp = interpolate(bgFullReshaped,(NoInterp(),BSpline(Linear(    )),NoInterp(),NoInterp(),NoInterp())); #Interpolire bgFrames    Full auf
     end
 else
   if zBG !=1
      itp = interpolate(bgFullReshaped,(NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),NoInterp())); #Interpolire bgFrames    Full auf
   else
      itp = interpolate(bgFullReshaped,(NoInterp(),NoInterp(),NoInterp(),NoInterp(),NoInterp())); #Interpolire bgFrames        Full auf
   end
 end

  bgFramesFullInterp = itp(collect(1:xBG),range(1,yBG,length=yFG),range(1,zBG,length=zFG),collect(1:size(bgFullReshaped,4)),collect(1:size(bgFullReshaped,5)))

  params[:acqNumFrames] = numFGFrames + prod(size(bgFramesFullInterp)[1:3])

  tt = [round.(Int,collect(range(1,xBG+xFG,length=xBG))).+(kk-1)*(xBG+xFG) for kk =1:yFG*zFG];
  idxBGFrames = vcat(tt...);
  idxAllFrames = collect(1:params[:acqNumFrames])
  idxAllFrames[idxBGFrames] .= 0
  idxFGFrames = idxAllFrames[idxAllFrames.!=0];

  params[:measFramePermutation] =  vcat(idxFGFrames,idxBGFrames)
  params[:measIsBGFrame] = vcat(zeros(numFGFrames),ones(length(idxBGFrames)))
  dataTemp_ = vcat(fgFrames,reshape(bgFramesFullInterp,:,size(bgFullReshaped)[4:5]...,1))  #size(params[:measData])(22621, 26929, 3, 1)
  params[:measData] = convert(Array{Complex{Float32},4},dataTemp_);  #size(params[:measData])(22621, 26929, 3, 1)
  saveasMDF(filenameOut, params)
end


function convertCustomSF(filenameOut::String, f::Array{BrukerFileMeas}, fBG::BrukerFile,nAverages::Int64,skipSwitchingFrames::Int64; nAveragesBG = nAverages,skipSwitchingFramesBG=skipSwitchingFrames)

  paramsArray =[];
  paramsBG = loadMetadata(fBG)

  params = loadMetadata(f[1])
  loadMeasParams(f[1], params, skipMeasData = true)
  loadCalibParams(f[1], params)

  for fi in f[1:end]

  params_ = loadMetadata(fi)
  loadMeasParams(fi, params_, skipMeasData = true)
  loadCalibParams(fi, params_)
  push!(paramsArray,params_);
  end

  for i =2:length(f)
  params[:acqOffsetField] = cat(params[:acqOffsetField],paramsArray[i][:acqOffsetField],dims=3)
  end

  params[:calibSize] = [length(union(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  numFGFrames = prod(params[:calibSize])

  println(params[:calibSize])

  params[:acqGradient] = reshape(params[:acqGradient][:,:,1,1],3,3,1,1)
  params[:dfPhase] =  reshape(params[:dfPhase][:,:,1],1,3,1)
  params[:calibDeltaSampleSize] = [0.0, 0.0, 0.0] #Todo
  #params[:time] = now()
  params[:calibMethod] = "BrukerCustom"
  #params[:version] = v"2.0.0"
  params[:dfStrength] = reshape(params[:dfStrength][:,:,1,1],1,3,1)
  params[:experimentIsCalibration] = true
  #params[:uuid] = uuid4()
  params[:measIsFourierTransformed] = true
  params[:measIsFastFrameAxis] = true
  params[:calibFovCenter] = [mean(extrema(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  params[:acqNumPeriodsPerFrame] = 1

  params[:measIsFramePermutation] = 1

  if rxHasTransferFunction(f[1]) &&!measIsTFCorrected(f[1])
    params[:transferFunction] = sampleTF(TransferFunction("/opt/data/TFs/TFHDF5/TF_3DMouse_NM0_07102019.h5"), f[1])
#    params[:transferFunction] = rxTransferFunction(f[1])
  end

  params[:calibFov] = [sum(abs.(extrema(params[:acqOffsetField][ll,1,:])))./abs.(params[:acqGradient][ll,ll]) for ll in collect(1:3)]
  params[:calibOrder] = "xyz"

  params[:acqOffsetField] = reshape([mean(extrema(params[:acqOffsetField][ll,1,:])) for ll in collect(1:3)],3,1,1)
  params[:acqNumAverages] = nAverages

println("Part1")
  # Daten Laden
#fgFrames =zeros(ComplexF64,prod(Int64,params[:calibSize]),26929,3,1)
@time global fgFrames =loadAndProcessFFData(f[1],nAverages, skipSwitchingFrames)
#global tmp =1
for (ii,ffi)  in enumerate(f[2:end])
println(ffi)
#pos = prod(Int64,[length(union(paramsArray[ii][:acqOffsetField][ll,1,:])) for ll in collect(1:3)])
#@time fgFrames[tmp:tmp+pos-1,:,:,:] = loadAndProcessFFData(ffi,nAverages, skipSwitchingFrames) # allePos x freqs x r x 1
@time global fgFrames = cat(fgFrames,loadAndProcessFFData(ffi,nAverages, skipSwitchingFrames),dims=1) # allePos x freqs x r x 1
#global tmp = tmp+pos
end

@time  bgFramesFull = loadAndProcessFFData(fBG,nAveragesBG, skipSwitchingFramesBG) # allePos x freqs x r x 1
  nAveragesBGHalf = nAveragesBG >1 ? div(nAveragesBG,2) : nAveragesBG
@time  bgFramesHalf = loadAndProcessFFData(fBG,div(nAveragesBG,2), skipSwitchingFramesBG+nAveragesBG-div(nAveragesBG,2)) # allePos x freqs x r x 1

  params[:calibSNR] = calculateSNRCustomSF(f[1],fgFrames,bgFramesFull,bgFramesHalf)

  (xBG,yBG,zBG) = [length(union(paramsBG[:acqOffsetField][ll,1,:])) for ll in collect(1:3)]
  (xFG,yFG,zFG) = params[:calibSize]

  bgFullReshaped = reshape(bgFramesFull,xBG,yBG,zBG,rxNumFrequencies(fBG),rxNumChannels(fBG))

  itp = interpolate(bgFullReshaped,(NoInterp(),BSpline(Linear()),BSpline(Linear()),NoInterp(),NoInterp())); #Interpolire bgFramesFull auf

  bgFramesFullInterp = itp(collect(1:xBG),range(1,yBG,length=yFG),range(1,zBG,length=zFG),collect(1:size(bgFullReshaped,4)),collect(1:size(bgFullReshaped,5)))

  params[:acqNumFrames] = numFGFrames + prod(size(bgFramesFullInterp)[1:3])

  tt = [round.(Int,collect(range(1,xBG+xFG,length=xBG))).+(kk-1)*(xBG+xFG) for kk =1:yFG*yFG];
  idxBGFrames = vcat(tt...);
  idxAllFrames = collect(1:params[:acqNumFrames])
  idxAllFrames[idxBGFrames] .= 0
  idxFGFrames = idxAllFrames[idxAllFrames.!=0];

  params[:measFramePermutation] =  vcat(idxFGFrames,idxBGFrames)
  params[:measIsBGFrame] = vcat(zeros(numFGFrames),ones(length(idxBGFrames)))
  dataTemp_ = vcat(fgFrames,reshape(bgFramesFullInterp,:,size(bgFullReshaped)[4:5]...,1))  #size(params[:measData])(22621, 26929, 3, 1)
  params[:measData] = convert(Array{Complex{Float32},4},dataTemp_);  #size(params[:measData])(22621, 26929, 3, 1)
  saveasMDF(filenameOut, params)
end

"""
This function returns the time data `timeData` Array(Float32,numPos,samplingPoints,channels,1,numFrames)
of an FF-measurement with `numPos` as number of FF-positions. The function splits the entire FF-measurement
in measurement per FF-position. As a second return argument the corresponding positions vector
`pos` Array(Float64,3,1,numPos) is returned.
"""
function getFFdataPerPos(f::MPIFile, nAverages::Int64, skipSwitchingFrames::Int64; addToEnd=0)
    dataFilename = joinpath(f.path,"rawdata.job0")
    AllFrames = acqNumPeriodsPerFrame(f)
    numPos_= size(acqOffsetField(f),3)
    numPos = convert(Int64,div(numPos_,nAverages+skipSwitchingFrames))
    timeData_ = zeros(Float32,numPos,rxNumSamplingPoints(f),rxNumChannels(f),1,nAverages)

    posFF = acqOffsetField(f).*[-1,-1,-1] #Field to position
    pos = zeros(Float64,3,1,numPos)

    measTime = measData(f)

    for p = collect(1:numPos)
        st = p*skipSwitchingFrames+(p-1)*(nAverages+addToEnd)+1
        ind = collect(st:st+nAverages-1)
        pos_ = unique(posFF[:,:,ind],dims=3)
        if size(pos_,3) != 1
            error("FF Positions not unqiue! Check Average and skipSwitchingFrames Parameter!")
        end
        pos[:,:,p] = pos_
        tempTimeData = measTime[:,:,ind,:]
        timeData_[p,:,:,:,:] = permutedims(tempTimeData, [1,2,4,3])# switch last two dimensions: frames and patches switch
    end
    timeData = convert(Array{Int16,5},timeData_)
    return timeData, pos
end

"""
This function spilts a FF-measurement according to it FF-positions in measurements per position
and treats them as single measurements. It returns `fileNames` Array{String,1}
and the correspodning `paramsArray` as Array{Dict,1}, which are the data dictionaries.
Both can be used to save the data with saveasMDF(...).
"""
function prepareAsMDFSingleMeasurement(f::MPIFile, numFrames::Int64, skipSwitchingFrames::Int64; addToEnd=0)
    global params = loadMetadata(f)

    # get FF data
    timeDataPerPos, pos = getFFdataPerPos(f, numFrames, skipSwitchingFrames, addToEnd=addToEnd);

    # adjust meta data
    params[:acqGradient] = acqGradient(f)[:,:,1:1,1:1]
    params[:acqNumPeriodsPerFrame] = 1
    params[:acqNumFrames] = numFrames
    params[:dfStrength] = dfStrength(f)[:,:,1:1]
    params[:dfPhase] = dfPhase(f)[:,:,1:1]

    # adjust data and meta data per position
    paramsArray = Array{Any,1}(undef,size(pos,3))
    fileNames = Array{String,1}(undef,size(pos,3))
    for p=1:size(pos,3)
        tempParas=deepcopy(params)
        filename = string(params[:experimentNumber],"_pos",p)
        tempParas[:acqOffsetField] = pos[:,:,p:p] # all
        tempParas[:measIsFrequencySelection] = false
        tempParas[:measIsBGCorrected] = false
        tempParas[:measIsFramePermutation] = false
        tempParas[:measIsFastFrameAxis] = false
        tempParas[:measIsFourierTransformed] = false
        tempParas[:measIsSpectralLeakageCorrected] = false
        tempParas[:measIsBGFrame] = convert(Array{Bool,1},zeros(size(timeDataPerPos,5)));
        tempParas[:measIsTFCorrected] = false
        tempParas[:measData] = timeDataPerPos[p,:,:,:,:];

        fileNames[p] = filename
        paramsArray[p] = tempParas
    end
    return fileNames, paramsArray
end

hasKeyAndValue(paramDict,param) = haskey(paramDict, param) && paramDict[param] != nothing

function writeIfAvailable(file, paramOut, paramDict, paramIn )
  if hasKeyAndValue(paramDict, paramIn)
    write(file, paramOut, paramDict[paramIn])
  end
end


function saveasMDF(file::HDF5.File, params::Dict{T,Any}) where  {T<:AbstractString}
  paramsSymbol = Dict{Symbol,Any}([Symbol(k)=>v for (k,v) in params])
  saveasMDF(file, paramsSymbol)
end

function saveasMDF(file::HDF5.File, params::Dict{Symbol,Any})
  # general parameters
  write(file, "/version", "2.0.1")
  write(file, "/uuid", string(get(params,:uuid,uuid4() )))
  write(file, "/time", "$( get(params,:time, Dates.unix2datetime(time())) )")

  # study parameters
  write(file, "/study/name", get(params,:studyName,"default") )
  write(file, "/study/number", get(params,:studyNumber,0))
  if hasKeyAndValue(params,:studyUuid)
    studyUuid = params[:studyUuid]
  else
    studyUuid = uuid4()
  end
  write(file, "/study/uuid", string(studyUuid))
  write(file, "/study/description", get(params,:studyDescription,"n.a."))
  if hasKeyAndValue(params,:studyTime)
    write(file, "/study/time", string(params[:studyTime]))
  end

  # experiment parameters
  write(file, "/experiment/name", get(params,:experimentName,"default") )
  write(file, "/experiment/number", get(params,:experimentNumber,0))
  if hasKeyAndValue(params,:experimentUuid)
    expUuid = params[:experimentUuid]
  else
    expUuid = uuid4()
  end
  write(file, "/experiment/uuid", string(expUuid))
  write(file, "/experiment/description", get(params,:experimentDescription,"n.a."))
  write(file, "/experiment/subject", get(params,:experimentSubject,"n.a."))
  write(file, "/experiment/isSimulation", Int8(get(params,:experimentIsSimulation,false)))

  # tracer parameters
  write(file, "/tracer/name", get(params,:tracerName,"n.a") )
  write(file, "/tracer/batch", get(params,:tracerBatch,"n.a") )
  write(file, "/tracer/vendor", get(params,:tracerVendor,"n.a") )
  write(file, "/tracer/volume", get(params,:tracerVolume,0.0))
  write(file, "/tracer/concentration", get(params,:tracerConcentration,0.0) )
  write(file, "/tracer/solute", get(params,:tracerSolute,"Fe") )
  tr = [string(t) for t in get(params,:tracerInjectionTime, [Dates.unix2datetime(time())]) ]
  write(file, "/tracer/injectionTime", tr)

  # scanner parameters
  write(file, "/scanner/facility", get(params,:scannerFacility,"n.a") )
  write(file, "/scanner/operator", get(params,:scannerOperator,"n.a") )
  write(file, "/scanner/manufacturer", get(params,:scannerManufacturer,"n.a"))
  write(file, "/scanner/name", get(params,:scannerName,"n.a"))
  write(file, "/scanner/topology", get(params,:scannerTopology,"FFP"))

  # acquisition parameters
  write(file, "/acquisition/numAverages",  params[:acqNumAverages])
  write(file, "/acquisition/numFrames", get(params,:acqNumFrames,1))
  write(file, "/acquisition/numPeriodsPerFrame", get(params,:acqNumPeriodsPerFrame,1))
  write(file, "/acquisition/startTime", "$( get(params,:acqStartTime, Dates.unix2datetime(time())) )")

  writeIfAvailable(file, "/acquisition/gradient", params, :acqGradient)
  writeIfAvailable(file, "/acquisition/offsetField", params, :acqOffsetField)

  # drivefield parameters
  write(file, "/acquisition/drivefield/numChannels", size(params[:dfStrength],2) )
  write(file, "/acquisition/drivefield/strength", params[:dfStrength])
  write(file, "/acquisition/drivefield/phase", params[:dfPhase])
  write(file, "/acquisition/drivefield/baseFrequency", params[:dfBaseFrequency])
  write(file, "/acquisition/drivefield/divider", params[:dfDivider])
  write(file, "/acquisition/drivefield/cycle", params[:dfCycle])
  if !haskey(params, :dfWaveform)
    params[:dfWaveform] = "sine"
  end
  write(file, "/acquisition/drivefield/waveform", params[:dfWaveform])

  # receiver parameters
  write(file, "/acquisition/receiver/numChannels", params[:rxNumChannels])
  write(file, "/acquisition/receiver/bandwidth", params[:rxBandwidth])
  write(file, "/acquisition/receiver/numSamplingPoints", params[:rxNumSamplingPoints])
  if !haskey(params, :rxUnit)
    params[:rxUnit] = "V"
  end
  write(file, "/acquisition/receiver/unit",  params[:rxUnit])
  write(file, "/acquisition/receiver/dataConversionFactor",  params[:rxDataConversionFactor])
  if hasKeyAndValue(params,:rxTransferFunction)
    tf = params[:rxTransferFunction]
    group = file["/acquisition/receiver"]
    writeComplexArray(group, "transferFunction", tf)
  end
  writeIfAvailable(file, "/acquisition/receiver/transferFunctionFileName", params, :rxTransferFunctionFileName)
  writeIfAvailable(file, "/acquisition/receiver/inductionFactor",  params, :rxInductionFactor)

  # measurements
  if hasKeyAndValue(params, :measData)
    meas = params[:measData]
    if eltype(meas) <: Complex
      group = create_group(file,"/measurement")
      writeComplexArray(group, "/measurement/data", meas)
    else
      write(file, "/measurement/data", meas)
    end
    write(file, "/measurement/isFourierTransformed", Int8(params[:measIsFourierTransformed]))
    write(file, "/measurement/isSpectralLeakageCorrected", Int8(params[:measIsSpectralLeakageCorrected]))
    write(file, "/measurement/isTransferFunctionCorrected", Int8(params[:measIsTFCorrected]))
    write(file, "/measurement/isFrequencySelection", Int8(params[:measIsFrequencySelection]))
    write(file, "/measurement/isBackgroundCorrected",  Int8(params[:measIsBGCorrected]))
    write(file, "/measurement/isFastFrameAxis",  Int8(params[:measIsFastFrameAxis]))
    write(file, "/measurement/isFramePermutation",  Int8(params[:measIsFramePermutation]))
    writeIfAvailable(file, "/measurement/frequencySelection",  params, :measFrequencySelection)

    if hasKeyAndValue(params, :measFramePermutation)
      write(file, "/measurement/framePermutation", params[:measFramePermutation] )
    end
    if hasKeyAndValue(params, :measIsBGFrame)
      write(file, "/measurement/isBackgroundFrame", convert(Array{Int8}, params[:measIsBGFrame]) )
    end
    if hasKeyAndValue(params, :measIsSparsityTransformed)
      write(file, "/measurement/isSparsityTransformed", params[:measIsSparsityTransformed] )
      write(file, "/measurement/subsamplingIndices", params[:measSubsamplingIndices] )
      write(file, "/measurement/sparsityTransformation", params[:measSparsityTransformation] )
    end
  end

  # calibrations
  writeIfAvailable(file, "/calibration/snr",  params, :calibSNR)
  writeIfAvailable(file, "/calibration/fieldOfView",  params, :calibFov)
  writeIfAvailable(file, "/calibration/fieldOfViewCenter",  params, :calibFovCenter)
  writeIfAvailable(file, "/calibration/size",  params, :calibSize)
  writeIfAvailable(file, "/calibration/order",  params, :calibOrder)
  writeIfAvailable(file, "/calibration/positions",  params, :calibPositions)
  writeIfAvailable(file, "/calibration/offsetField",  params, :calibOffsetField)
  writeIfAvailable(file, "/calibration/deltaSampleSize",  params, :calibDeltaSampleSize)
  writeIfAvailable(file, "/calibration/method",  params, :calibMethod)
  if hasKeyAndValue(params, :calibIsMeanderingGrid)
    write(file, "/calibration/isMeanderingGrid", Int8(params[:calibIsMeanderingGrid]))
  end
  writeIfAvailable(file, "/calibration/_temperatures",  params, :calibTemperatures)
  # reconstruction
  if hasKeyAndValue(params, :recoData)
    write(file, "/reconstruction/data", params[:recoData])
    write(file, "/reconstruction/fieldOfView", params[:recoFov])
    write(file, "/reconstruction/fieldOfViewCenter", params[:recoFovCenter])
    write(file, "/reconstruction/size", params[:recoSize])
    write(file, "/reconstruction/order", get(params, :recoOrder, "xyz"))
    if hasKeyAndValue(params, :recoPositions)
      write(file, "/reconstruction/positions", params[:recoPositions])
    end
    if hasKeyAndValue(params, :recoParameters)
      saveParams(file, "/reconstruction/_parameters", params[:recoParameters])
    end
  end

  writeIfAvailable(file, "/custom/auxiliaryData", params, :auxiliaryData)
end
