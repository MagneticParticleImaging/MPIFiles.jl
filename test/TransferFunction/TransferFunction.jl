using MPIFiles
using Test

@testset "Testing TransferFunction submodule" begin

dir = @__DIR__

a = TransferFunction(dir*"/example.s1p", frequencyWeighting=true)
if isfile(dir*"/example.h5")
  rm(dir*"/example.h5")
end
MPIFiles.save(dir*"/example.h5", a)
b = TransferFunction(dir*"/example.h5")

@test a.freq == b.freq
@test a.data == b.data
@test a[[1],1] == a[[0.0],1]

c = MPIFiles.combine(MPIFiles.combine(a,a),a)
rm(dir*"/example.h5")
MPIFiles.save(dir*"/example.h5", c)
@test a[[1],1] == c[[1],2]
@test a[[1],1] == c[[1],3]



fnMeasBruker = dir*"/../measurement_Bruker"
measBruker = MPIFile(fnMeasBruker)
tf = sampleTF(c, measBruker)
@test size(tf) == (rxNumFrequencies(measBruker), rxNumChannels(measBruker))

end
