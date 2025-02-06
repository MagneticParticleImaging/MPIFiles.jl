using MPIFiles, Test

@testset "Testing DatasetStore submodule" begin

# Bruker store

storeA = BrukerDatasetStore(joinpath(datadir, "BrukerStore"))
@test readonly(storeA) == true

studiesA = getStudies(storeA)
@test [s.name for s in studiesA] == ["SF", "Wuerfelphantom", "Size Distributions", "EasyAxisContrast"]

# MDF store

storeB = MDFDatasetStore(joinpath(tmpdir, "MDFStoreA"))
@test readonly(storeB) == false
empty!(storeB)

# First add studies and remove them
studyA = Study(storeB,"StudyA")
studyB = Study(storeB,"StudyB")

addStudy(storeB, studyA)
addStudy(storeB, studyB)

studies = getStudies(storeB)
@test studies == [studyA, studyB]

remove(studyB)
studies = getStudies(storeB)
@test studies == [studyA]

empty!(storeB)


# Now lets copy over the Bruker studies
exportData(storeA, storeB, keepExpNum=true)
storeC = MDFDatasetStore(joinpath(tmpdir, "MDFStoreB"))
empty!(storeC)
exportData(storeB, storeC, SNRThresh=4.0)
createArtifact(storeB, "https://")

@test validate(storeA)
@test validate(storeB)
@test validate(storeC)

studiesB = getStudies(storeB)
studiesC = getStudies(storeC)
@test getExperiments(studiesB[2])[1].num == 18 broken = true
@test getExperiments(studiesC[2])[1].num == 1 broken = true

# export into new/existing Study
#newStudy = Study(storeB, "NewStudy")
#exportData(studiesB[2], storeB, newStudy)
#@info [s.name for s in getStudies(storeB)] 
#@test getExperiments(newStudy)[1].num == 1


# Next we test changing parameters and study names.
# TODO broken
#s = studiesC[2]
#@test validate(s)
#changeStudy(s, "NewStudyName")
#studiesD = getStudies(storeC)
#sNew = studiesD[end]
#@test validate(sNew)
#@test sNew.name == "NewStudyName"
#enforceStudy(sNew)
#@test validate(sNew)


# Experiment handling

#exps = getExperiments(storeB, studyA)
#@test exps == Experiment[]
#@info getNewExperimentNum(storeB,studyA) == 1


# Calibration handling
@test getNewCalibNum(storeB) == 3

# Distributed Storage
@testset "Distributd Store" begin
  localStore = MDFDatasetStore(joinpath(tmpdir, "MDFStoreLocal"))
  empty!(localStore)
  remoteStore = DDatasetStore(joinpath(tmpdir, "MDFStoreLocal"); worker = 1)
  exportData(storeA, localStore, keepExpNum=true)

  @test readonly(localStore) == readonly(remoteStore)
  @test studydir(localStore) == studydir(remoteStore)
  @test calibdir(localStore) == calibdir(remoteStore)
  @test getCalibStudy(localStore).uuid == getCalibStudy(remoteStore).uuid
  @test getnewCalibNum(localStore) == getnewCalibNum(remoteStore)
  @test getNewCalibPath(localStore)[1:end-5] == getNewCalibPath(remoteStore)[1:end-5] # this command creates files so we skip "X.mdf"

  @test length(getStudies(localStore)) == length(getStudies(remoteStore))
  localStudy = getStudy(localStore, "Wuerfelphantom")
  remoteStudy = getStudy(remoteStore, "Wuerfelphantom")
  @test localStudy.uuid == remoteStudy.uuid
  @test path(localStudy) == path(remoteStudy)
  @test getNewExperimentNum(localStudy) == getNewExperimentNum(remoteStudy)
  @test getNewExperimentPath(localStudy)[1:end-5] == getNewExperimentPath(remoteStudy)[1:end-5] # this command creates files so we skip "X.mdf"

  # Adding Study
  newStudy = Study(remoteStore, "RemoteStudy")
  @test first(getStudies(localStore, "RemoteStudy")).uuid == newStudy.uuid
  remove(newStudy)
  @test isempty(getStudies(localStore, "RemoteStudy"))

  # Experiement
  localStudy = first(getStudies(localStore, "Wuerfelphantom"))
  remoteStudy = first(getStudies(remoteStore, "Wuerfelphantom"))
  @test getMDFStudyFolderName(remoteStudy) == getMDFStudyFolderName(localStudy)
  localExp = getExperiment(localStudy, 18)
  remoteExp = getExperiment(remoteStudy, 18)
  @test MPIFile(remoteExp) isa DMPIFile
  @test path(localExp) == path(remoteExp)
  @test localExp.time == remoteExp.time
  @test isnothing(getExperiment(remoteStudy, typemax(Int64)))
  @test length(getExperiments(remoteStudy)) == length(getExperiments(localStudy))
  @test MPIFiles.iscalib(remoteExp) == MPIFiles.iscalib(localExp)

  #empty!(remoteStore)
  #@test isempty(getStudies(localStore))
end


end