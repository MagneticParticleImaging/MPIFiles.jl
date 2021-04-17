using MPIFiles, Test

@testset "Testing DatasetStore submodule" begin

store = MDFDatasetStore(joinpath(@__DIR__, "MDFStore"))
@test readonly(store) == false
empty!(store)

# Study handling

studyA = Study(store,"StudyA")
studyB = Study(store,"StudyB")

addStudy(store, studyA)
addStudy(store, studyB)

studies = getStudies(store)
@test studies == [studyA, studyB]

remove(studyB)
studies = getStudies(store)
@test studies == [studyA]

# Experiment handling

exps = getExperiments(store, studyA)
@test exps == Experiment[]
@info getNewExperimentNum(store,studyA) == 1


# Calibration handling
@info getNewCalibNum(store) == 1

end