fnRecoV1 = "reco_V1.mdf"
fnRecoV2 = "reco_V2.mdf"
fnRecoV2b = "reco_V2b.mdf"

if !isfile(fnRecoV1)
  streamReco = get("http://media.tuhh.de/ibi/mdf/reco_V1.mdf")
  save(streamReco, fnRecoV1)
end

mdfv1 = MPIFile(fnRecoV1)
@test typeof(mdfv1) == MDFFileV1

saveasMDF(fnRecoV2, fnRecoV1)

mdfv2 = MPIFile(fnRecoV2)
@test typeof(mdfv2) == MDFFileV2

c1 = loadRecoDataMDF(fnRecoV1)
c2 = loadRecoDataMDF(fnRecoV2)
saveRecoDataMDF(fnRecoV2b,c1)

mdfv2b = MPIFile(fnRecoV2b)
@test typeof(mdfv2b) == MDFFileV2

for mdf in (mdfv1,mdfv2,mdfv2b)
  println("Test $mdf")
  @test recoSize(mdf) == [40,40,1]
  @test norm(recoFovCenter(mdf) - [0.0,0.0,0.0]) < eps()
  @test recoFov(mdf) == [0.044,0.044,0.001]
  @test recoOrder(mdf) == "xyz"
  @test recoPositions(mdf) == nothing

  @test size( recoData(mdf) ) == (1, 1600, 1)
end

@test c1 == c2
