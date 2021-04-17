using MPIFiles, Test

@testset "Testing DatasetStore submodule" begin

# Bruker store

storeA = BrukerDatasetStore(joinpath(@__DIR__, "data", "BrukerStore"))
@test readonly(storeA) == true

studiesA = getStudies(storeA)
@test id.(studiesA) == ["SF", "Wuerfelphantom", "Size Distributions", "EasyAxisContrast"]

# MDF store

storeB = MDFDatasetStore(joinpath(@__DIR__, "data", "MDFStoreA"))
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
storeC = MDFDatasetStore(joinpath(@__DIR__, "data", "MDFStoreB"))
empty!(storeC)
exportData(storeB, storeC)

# Experiment handling

#exps = getExperiments(storeB, studyA)
#@test exps == Experiment[]
#@info getNewExperimentNum(storeB,studyA) == 1


# Calibration handling
@info getNewCalibNum(storeB) == 1

end