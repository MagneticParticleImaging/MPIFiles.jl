using MPIFiles
using Test

fnMeasBruker = "./data/measurement"
tfpath = "./data/transferFunction/example.s1p"
tfh5path = "./data/transferFunction/example.h5"

@testset "Testing TransferFunction submodule" begin
  a = TransferFunction(tfpath, frequencyWeighting=true)
  if isfile(tfh5path)
    rm(tfh5path)
  end
  MPIFiles.save(tfh5path, a)
  b = TransferFunction(tfh5path)

  @test a.freq == b.freq
  @test a.data == b.data
  @test a[[1],1] == a[[0.0],1]

  c = MPIFiles.combine(MPIFiles.combine(a,a),a)
  rm(tfh5path)
  MPIFiles.save(tfh5path, c)
  @test a[[1],1] == c[[1],2]
  @test a[[1],1] == c[[1],3]

  measBruker = MPIFile(fnMeasBruker)
  tf = sampleTF(c, measBruker)
  @test size(tf) == (rxNumFrequencies(measBruker), rxNumChannels(measBruker))
end
