@testset "Testing Cartesian Sequence submodule" begin

fnSMTD = "./data/mdf/systemMatrixCartesian.mdf"

smTD = MPIFile(fnSMTD)
@test typeof(smTD) <: MDFFileV2


@test acqNumPeriodsPerFrame(smTD) == 6500
@test size(getSystemMatrix(smTD,1:10)) == (81,10,6500)
@test size(getMeasurements(smTD)) == (76, 2, 6500, 81)
@test size(getMeasurements(smTD, numPeriodAverages=65)) == (76, 2, 100, 81)
@test size(getMeasurements(smTD, numPeriodAverages=65, numPeriodGrouping=100)) == (7600, 2, 1, 81)
@test size(getMeasurementsFD(smTD, frequencies=1:10)) == (10, 6500, 81)

end