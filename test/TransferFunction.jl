using MPIFiles
using Test

fnMeasBruker = joinpath(datadir,"BrukerStore","20150915_102110_Wuerfelphantom_1_1","18")
fnMeasMDFv1 = joinpath(datadir, "mdf", "measurement_V1.mdf")
fnTmpMDFv1 = joinpath(tmpdir, "mdf", "measurement_V1_tmp.mdf")
fnMeasMDFv2 = joinpath(tmpdir,"mdfim","measurement_V2.mdf")
tfpath = joinpath(datadir,"transferFunction","example.s1p")
tfh5pathTmp = joinpath(tmpdir,"transferFunction","example.h5")
tfh5path = joinpath(datadir,"transferFunction","tf.h5")

@testset "Testing TransferFunction submodule" begin
  a = TransferFunction(tfpath, frequencyWeighting=true)
  if isfile(tfh5pathTmp)
    rm(tfh5pathTmp)
  end
  MPIFiles.save(tfh5pathTmp, a)
  b = TransferFunction(tfh5pathTmp)

  @test a.freq == b.freq
  @test a.data == b.data
  @test a[1,1] == a(0.0,1)

  c = MPIFiles.combine(MPIFiles.combine(a,a),a)
  rm(tfh5pathTmp)
  MPIFiles.save(tfh5pathTmp, c)
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

  tf = TransferFunction(f, data)
  tf_kHz = TransferFunction(f*u"kHz", data)
  @test tf(100) ≈ tf_kHz(100e3)

  data_units = [1u"V/T" ./(1 .+im*f/1e4) 1u"V/(A*m^2)" ./(1 .+im*f/1e3)]
  tf = TransferFunction(f, data_units)
  @test_throws ErrorException TransferFunction([0,1], [1.0u"V/T" 1.0u"V/A"; 1.0u"V/T" 1.0u"m"])
  
  tf_comp = TransferFunction([0,1,2], [1,im*2,-im*3])
  tf_deg = TransferFunction([0,1,2], [1,2,3], [0,90,270]u"°")
  tf_rad = TransferFunction([0,1,2], [1,2,3], [0,π/2,-π/2])
  @test all(tf_comp([0,1,2]) .≈ tf_deg([0,1,2]))
  @test all(tf_deg([0,1,2]) .≈ tf_rad([0,1,2]))

  f1 = collect(range(0,1e6,step=1e3));
  data1 = 1 ./ (1 .+ im*f1/1e4 );

  f2 = collect(range(0,1e7,step=1e4));
  data2 = 1 ./ (1 .+ im*f2/1e4 );

  tf1 = TransferFunction(f1, data1)
  tf2 = TransferFunction(f2, data2)

  @test_throws ErrorException tf_c = combine(tf1,tf2)
  tf_c = combine(tf1,tf2,interpolate=true)
  @test tf_c(101.5e3, :) ≈ [tf1(101.5e3) tf2(101.5e3)] atol=1e-12

  # prepare files
  cp(fnMeasMDFv1, fnTmpMDFv1, force=true)
  chmod(fnTmpMDFv1, 0o777)
  f_mdfv1 = MPIFile(fnMeasMDFv1)
  f_mdfv1_tmp = MPIFile(fnTmpMDFv1)
  f_mdfv2 = MPIFile(fnMeasMDFv2)
  
  
  # create TF from file
  tf = TransferFunction(f_mdfv2)
  @test sampleTF(tf, f_mdfv2) ≈ rxTransferFunction(f_mdfv2)

  # set TF to file that has no TF
  @test_throws ErrorException getMeasurementsFD(f_mdfv1, tfCorrection=true)
  newtf = TransferFunction([0,10e6],[1.0+0.0*im 1.0+0.0*im 1.0+0.0*im; 1.0+0.0*im 1.0+0.0*im 1.0+0.0*im])
  setTF(f_mdfv1_tmp, newtf)
  @test rxHasTransferFunction(f_mdfv1_tmp)
  @test getMeasurementsFD(f_mdfv1_tmp, tfCorrection=true) == getMeasurementsFD(f_mdfv1_tmp, tfCorrection=false)

  # read TF from h5 file
  tf = TransferFunction(tfh5path)
  # set TF using file path
  setTF(f_mdfv1_tmp, tfh5path)
  @test rxTransferFunction(f_mdfv1_tmp) == sampleTF(tf, f_mdfv1_tmp)
  # read TF from h5 file with channel selection
  tf2 = TransferFunction(tfh5path, channels=[1,3])
  tf3 = TransferFunction(tfh5path, channels=2)
  @test tf[:, 1] == tf2[:, 1]
  @test tf[:, 3] == tf2[:, 2]
  @test tf[:, 2] == tf3[:, 1]

  setTF(MDFv2InMemory(f_mdfv2), newtf)

  # test printing
  out = sprint() do io
    show(io, MIME"text/plain"(),newtf)
  end
  @test out == "MPIFiles.TransferFunction: \n\t3 channel(s), units of [\"NoUnits\", \"NoUnits\", \"NoUnits\"]\n\t2 frequency samples from 0.0 Hz to 1.0e7 Hz"


  # loading and processing tf
  tf1 = MPIFiles.load_tf_fromVNA(tfpath, R=50, N=8, A=1e-3^2*pi)
  tf2 = MPIFiles.load_tf_fromVNA(tfpath, R=50, N=8, d=2e-3)
  tf3 = MPIFiles.load_tf_fromVNA(tfpath, R=50, N=8, r=1e-3)
  @test tf1.data == tf2.data == tf3.data
  @test_throws ErrorException MPIFiles.load_tf_fromVNA(tfpath, R=50, d=2e-3)
  @test_throws ErrorException MPIFiles.load_tf_fromVNA(tfpath, R=50, N=8, d=2e-3, r=2e-3)
  MPIFiles.load_tf_fromVNA(tfpath, R=50, N=8, d=2e-3, frequencyWeighting=true)
end
