using MPIFiles, Test

@testset "Testing DatasetStore submodule" begin

# Bruker store

storeA = BrukerDatasetStore(joinpath(datadir, "BrukerStore"))
@test readonly(storeA) == true

studiesA = getStudies(storeA)
@test id.(studiesA) == ["SF", "Wuerfelphantom", "Size Distributions", "EasyAxisContrast"]

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
exportData(storeA, storeB)
storeC = MDFDatasetStore(joinpath(tmpdir, "MDFStoreB"))
empty!(storeC)
exportData(storeB, storeC, SNRThresh=4.0)
makeTarGzip(storeB)

# Experiment handling

#exps = getExperiments(storeB, studyA)
#@test exps == Experiment[]
#@info getNewExperimentNum(storeB,studyA) == 1


# Calibration handling
@info getNewCalibNum(storeB) 

end