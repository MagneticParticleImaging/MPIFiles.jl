@testset "Testing Reco submodule" begin

fnRecoV1 = "./data/mdf/reco_V1.mdf"
fnRecoV2 = "./data/mdf/reco_V2.mdf"
fnRecoV2b = "./data/mdf/reco_V2b.mdf"

mdfv1 = MPIFile(fnRecoV1)
@test typeof(mdfv1) <: MDFFileV1

saveasMDF(fnRecoV2, fnRecoV1)

mdfv2 = MPIFile(fnRecoV2)
@test typeof(mdfv2) <: MDFFileV2

c1 = loadRecoData(fnRecoV1)
c2 = loadRecoData(fnRecoV2)
saveRecoData(fnRecoV2b,c1)

mdfv2b = MPIFile(fnRecoV2b)
@test typeof(mdfv2b) <: MDFFileV2

for mdf in (mdfv1,mdfv2,mdfv2b)
  @info "Test $mdf"
  @test recoSize(mdf) == [40,40,1]
  @test norm(recoFovCenter(mdf) - [0.0,0.0,0.0]) < eps() # Image conversion fails
  @test recoFov(mdf) == [0.044,0.044,0.001]
  @test recoOrder(mdf) == "xyz"
  @test recoPositions(mdf) == nothing

  @test size( recoData(mdf) ) == (1, 1600, 1)
end

@test c1 == c2

end
