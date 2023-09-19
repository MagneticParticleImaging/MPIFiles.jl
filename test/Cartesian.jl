@testset "Testing Cartesian Sequence submodule" begin

fnSMTD = joinpath(datadir,"mdf","systemMatrixCartesian.mdf")
fnSMFDMP = joinpath(tmpdir,"mdf","systemMatrixCartesianMP.mdf")
fnSMFDSP = joinpath(tmpdir,"mdf","systemMatrixCartesianSP.mdf")

smTD = MPIFile(fnSMTD)
@test typeof(smTD) <: MDFFileV2

freqs = collect(vec(CartesianIndices((10, 1))))
@test acqNumPeriodsPerFrame(smTD) == 6500
@test rxNumSamplingPoints(smTD) == 76
@test size(getSystemMatrix(smTD,freqs)) == (81,6500*10)
@test size(getMeasurements(smTD)) == (76, 2, 6500, 81)
@test size(getMeasurements(smTD, numPeriodAverages=65)) == (76, 2, 100, 81)
@test size(getMeasurements(smTD, numPeriodAverages=65, numPeriodGrouping=100)) == (7600, 2, 1, 81)
@test size(getMeasurementsFD(smTD, frequencies=freqs)) == (10, 6500, 81)

saveasMDF(fnSMFDMP, smTD, numPeriodAverages=65, applyCalibPostprocessing=true)
smFDMP = MPIFile(fnSMFDMP)

@test acqNumPeriodsPerFrame(smFDMP) == 100
@test rxNumSamplingPoints(smFDMP) == 76
@test size(getSystemMatrix(smFDMP,freqs)) == (81,10*100)

saveasMDF(fnSMFDSP, smTD, numPeriodAverages=65, applyCalibPostprocessing=true, numPeriodGrouping=100)
smFDSP = MPIFile(fnSMFDSP)

@test acqNumPeriodsPerFrame(smFDSP) == 1
@test rxNumSamplingPoints(smFDSP) == 7600
@test size(getSystemMatrix(smFDSP,freqs)) == (81,10)

end