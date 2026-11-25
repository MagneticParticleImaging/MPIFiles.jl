@testset "Recovery from broken file handles" begin
  # Simulates a broken handle (e.g. after an interrupted network connection) by closing the file descriptor behind HDF5's back
  function breakHandle!(file::HDF5.File)
    fapl = HDF5.API.h5f_get_access_plist(file)
    fd = try
      unsafe_load(convert(Ptr{Cint}, HDF5.API.h5f_get_vfd_handle(file, fapl)))
    finally
      HDF5.API.h5p_close(fapl)
    end
    # HDF5_jll on Windows uses the file descriptors of msvcrt
    ret = Sys.iswindows() ? ccall((:_close, "msvcrt"), Cint, (Cint,), fd) : ccall(:close, Cint, (Cint,), fd)
    ret == 0 || error("Could not close file descriptor $fd")
  end

  mkpath(joinpath(tmpdir, "brokenHandles"))
  fnV1 = joinpath(tmpdir, "brokenHandles", "measurement_V1.mdf")
  cp(joinpath(datadir, "mdf", "measurement_V1.mdf"), fnV1, force = true)
  fnV2 = joinpath(tmpdir, "brokenHandles", "measurement_V2.mdf")
  saveasMDF(fnV2, defaultMDFv2InMemory())
  versV2 = MPIFile(version, fnV2)

  @testset "Single handle" begin
    for fn in (fnV1, fnV2)
      f = MPIFile(fn)
      vers = version(f)
      data = isnothing(f.mmap_measData) ? nothing : copy(f.mmap_measData)
      oldFile = f.file
      breakHandle!(f.file)
      @test version(f) == vers
      @test f.file !== oldFile
      @test isnothing(data) || f.mmap_measData == data
      @test experimentHasMeasurement(f) isa Bool
      close(f)
      @test isnothing(f.mmap_measData)
    end
  end

  @testset "Multiple MDFFiles" begin
    f1 = MPIFile(fnV2)
    f2 = MDFFile(fnV2)
    breakHandle!(f1.file)
    @test version(f1) == versV2
    @test version(f2) == versV2
    f3 = MPIFile(fnV2)
    @test version(f3) == versV2
    close.((f1, f2, f3))
  end

  @testset "Opening a file with a broken handle" begin
    f1 = MPIFile(fnV2)
    breakHandle!(f1.file)
    f2 = MPIFile(fnV2)
    @test version(f2) == versV2
    @test version(f1) == versV2
    breakHandle!(f1.file)
    f3 = MDFFile(fnV2)
    @test version(f3) == versV2
    @test version(f1) == versV2
    close.((f1, f2, f3))
  end

  @testset "Handles outside of MDFFiles" begin
    f1 = MPIFile(fnV2)
    h = h5open(fnV2, "r")
    breakHandle!(f1.file)
    # HDF5 returns wrong data from its caches after the first error, therefore the file stays inaccessible
    for _ in 1:2
      @test_throws "close all other handles of this file" version(f1)
    end
    @test_throws "close all other handles of this file" MPIFile(fnV2)
    @test_throws "close all other handles of this file" MDFFile(fnV2)
    close(h)
    @test version(f1) == versV2
    f2 = MPIFile(fnV2)
    @test version(f2) == versV2
    close.((f1, f2))
  end

  @testset "Closing a broken handle" begin
    f = MPIFile(fnV2)
    breakHandle!(f.file)
    @test_throws HDF5.API.H5Error close(f)
    @test !isopen(f.file)
    GC.gc() # The finalizer of the broken handle must not try to close it again
    f = MPIFile(fnV2)
    @test version(f) == versV2
    close(f)
  end

  @testset "Manual reopen" begin
    f = MPIFile(fnV2)
    @test reopen!(f) === f
    @test version(f) == versV2
    close(f)
    @test reopen!(f) === f
    @test version(f) == versV2
    close(f)
  end
end
