@testset "Recovery from broken file handles" begin
  # Simulates a broken handle (e.g. after an interrupted network connection) by closing the file descriptor behind HDF5's back
  function breakHandle!(file::HDF5.File)
    fapl = HDF5.API.h5f_get_access_plist(file)
    fd = try
      unsafe_load(convert(Ptr{Cint}, HDF5.API.h5f_get_vfd_handle(file, fapl)))
    finally
      HDF5.API.h5p_close(fapl)
    end
    ccall(:close, Cint, (Cint,), fd)
  end

  mkpath(joinpath(tmpdir, "brokenHandles"))
  fnV1 = joinpath(tmpdir, "brokenHandles", "measurement_V1.mdf")
  cp(joinpath(datadir, "mdf", "measurement_V1.mdf"), fnV1, force = true)
  fnV2 = joinpath(tmpdir, "brokenHandles", "measurement_V2.mdf")
  saveasMDF(fnV2, defaultMDFv2InMemory())

  if Sys.iswindows()
    @test_skip "Simulating broken handles is not supported on Windows"
  else
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
      end
    end

    @testset "Multiple handles" begin
      f1 = MPIFile(fnV2)
      f2 = MPIFile(fnV2)
      vers = version(f1)
      breakHandle!(f1.file)
      @test_throws r"close\(\) all other MPIFiles" version(f1)
      @test_throws r"close\(\) all other MPIFiles" version(f2)
      @test_throws r"close\(\) all other MPIFiles" MPIFile(fnV2)
      @test_throws r"close\(\) all other MPIFiles" MDFFile(fnV2)
      close(f2)
      @test version(f1) == vers
      f3 = MPIFile(fnV2)
      @test version(f3) == vers
      close(f3)
      close(f1)
    end

    @testset "Manual reopen" begin
      f = MPIFile(fnV2)
      vers = version(f)
      @test reopen!(f) === f
      @test version(f) == vers
      close(f)
    end
  end
end
