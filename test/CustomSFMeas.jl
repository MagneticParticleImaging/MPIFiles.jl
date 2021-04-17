@testset "Testing Custom SF and Meas submodule" begin

customSF = "./data/BrukerStore/20200603_172811_EasyAxisContrast_1_1/245"
customMeas = "./data/BrukerStore/20200603_172811_EasyAxisContrast_1_1/245"

fnCustomSF_FF = "./data/mdf/customSystemMatrixFF.mdf"

fSF = MPIFile(customSF)
fMeas = MPIFile(customMeas)

data = MPIFiles.loadAndProcessFFData(fSF, 100, 4)
@test size(data) == (9,817,3,1)

convertCustomSF(fnCustomSF_FF, fSF, fSF, 100, 4)
dataSet = loadDataset(MPIFile(fnCustomSF_FF))
@test size(dataSet[:measData]) == (18,817,3,1)

timeData, pos = getFFdataPerPos(fMeas, 100, 4)
@test size(timeData) == (3,1632,3,1,100)
@test size(pos) == (3,1,3)

fileNames, paramsArray = prepareAsMDFSingleMeasurement(fMeas, 100, 4)
@test length(fileNames) == 3
@test length(paramsArray) == 3
end