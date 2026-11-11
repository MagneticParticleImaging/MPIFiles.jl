@testset "Forward reference from Bruker directory to MDF" begin
  # A Bruker experiment directory can contain a file `mdf` with the path
  # of an MDF file that should be opened instead of the Bruker data.
  mdfFile = joinpath(tmpdir, "mdf", "forwardRefTarget.mdf")
  saveasMDF(mdfFile, defaultMDFv2InMemory())

  @testset "Existing MDF file" begin
    brukerDir = mkpath(joinpath(tmpdir, "forwardRef", "valid"))
    write(joinpath(brukerDir, "mdf"), mdfFile)
    f = MPIFile(brukerDir)
    @test f isa MDFFileV2
    @test filepath(f) == mdfFile
    close(f)
  end

  @testset "Missing MDF file" begin
    brukerDir = mkpath(joinpath(tmpdir, "forwardRef", "broken"))
    write(joinpath(brukerDir, "mdf"), joinpath(tmpdir, "forwardRef", "doesNotExist.mdf"))
    f = @test_logs (:warn, r"forward ref") MPIFile(brukerDir)
    @test f isa BrukerFile
    @test f.path == brukerDir
  end
end
