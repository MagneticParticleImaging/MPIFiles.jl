using MPIFiles
using Test

fnMeasBruker = joinpath(datadir,"BrukerStore","20150915_102110_Wuerfelphantom_1_1","18")
tfpath = joinpath(datadir,"transferFunction","example.s1p")
tfh5path = joinpath(tmpdir,"transferFunction","example.h5")

@testset "Testing TransferFunction submodule" begin
  a = TransferFunction(tfpath, frequencyWeighting=true)
  if isfile(tfh5path)
    rm(tfh5path)
  end
  MPIFiles.save(tfh5path, a)
  b = TransferFunction(tfh5path)

  @test a.freq == b.freq
  @test a.data == b.data
  @test a[1,1] == a(0.0,1)

  c = MPIFiles.combine(MPIFiles.combine(a,a),a)
  rm(tfh5path)
  MPIFiles.save(tfh5path, c)
  @test a[1,1] == c[1,2]
  @test a[1,1] == c[1,3]

  measBruker = MPIFile(fnMeasBruker)
  tf = sampleTF(c, measBruker)
  @test size(tf) == (rxNumFrequencies(measBruker), rxNumChannels(measBruker))

  f = collect(range(0,1e6,step=1e3));
  data = 1 ./ (1 .+ im*f/1e4 );
  tf = TransferFunction(f, data) # TransferFunction of a simple lowpass filter
  @test tf[1] == 1.0
  @test tf[1:6,1] == data[1:6]
  @test tf(0) == tf(0,1)
  
  data = [1 ./(1 .+im*f/1e4) 1 ./(1 .+im*f/1e3)]
  tf = TransferFunction(f, data)
  @test tf[11,1] == data[11,1]
  @test tf[11,2] == data[11,2]
  @test tf[11,:] == data[11,:]
  @test tf(1e4, [1,2]) == [data[11,1] 1 ./(1 .+im*1e4/1e3)]
  @test tf(1e4,:) == tf(1e4, [1,2])
  @test tf(1e4,:) == tf(1e4, 1:2)
  @test tf(0) == tf(0,1)

  tf = TransferFunction(f, data, units=["V/V", "A/V"])
  @test tf(0) == 1u"V/V"
  @test tf(0,2) == 1u"A/V"


  f1 = collect(range(0,1e6,step=1e3));
  data1 = 1 ./ (1 .+ im*f1/1e4 );

  f2 = collect(range(0,1e7,step=1e4));
  data2 = 1 ./ (1 .+ im*f2/1e4 );

  tf1 = TransferFunction(f1, data1)
  tf2 = TransferFunction(f2, data2)

  @test_throws ErrorException tf_c = combine(tf1,tf2)
  tf_c = combine(tf1,tf2,interpolate=true)
  @test tf_c(101.5e3, :) ≈ [tf1(101.5e3) tf2(101.5e3)] atol=1e-12

end
