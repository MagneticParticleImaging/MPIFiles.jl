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

studiesB = getStudies(storeB)
studiesC = getStudies(storeC)
@test getExperiments(studiesB[2])[1].num == 18
@test getExperiments(studiesC[2])[1].num == 1

# Next we test changing parameters and study names.
s = studiesC[2]
@test validate(s)
changeStudy(s, "NewStudyName")
studiesD = getStudies(storeC)
sNew = studiesD[end]
@test validate(sNew)
@test sNew.name == "NewStudyName"
enforceStudy(sNew)
@test validate(sNew)


# Experiment handling

#exps = getExperiments(storeB, studyA)
#@test exps == Experiment[]
#@info getNewExperimentNum(storeB,studyA) == 1


# Calibration handling
@test getNewCalibNum(storeB) == 3

end