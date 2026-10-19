@testset "Project field" begin
  fnMeasBruker = joinpath(datadir,"BrukerStore","20150915_102110_Wuerfelphantom_1_1","18")
  testpath = joinpath(tmpdir, "mdf_project")
  measBruker = MPIFile(fnMeasBruker)
  
  fnMeasV2 = joinpath(testpath, "v2.mdf")
  saveasMDF(fnMeasV2, measBruker)

  @testset "v2" begin
    mdfv2 = MPIFile(fnMeasV2)

    # Initially we dont have any project information
    @test isnothing(projectName(mdfv2))
    @test isnothing(projectUuid(mdfv2))
    @test isnothing(projectTime(mdfv2))
    @test isnothing(projectDescription(mdfv2))

    # We can set these for in memory MDFs
    memV2 = MDFv2InMemory(mdfv2)
    testTitle = "Test"
    testDescription = "This is a Test"
    testDate = now()
    testUuid = uuid7()
    projectName(memV2, testTitle)
    @test projectName(memV2) == testTitle
    projectDescription(memV2, testDescription)
    @test projectDescription(memV2) == testDescription
    projectTime(memV2, testDate)
    @test projectTime(memV2) == testDate
    projectUuid(memV2, testUuid)
    @test projectUuid(memV2) == testUuid
    # We can persist these changes
    fnMeasMem = joinpath(testpath, "mem.mdf")
    saveasMDF(fnMeasMem, memV2)
    mdfFromMem = MPIFile(fnMeasMem)
    @test projectName(mdfFromMem) == testTitle
    @test projectDescription(mdfFromMem) == testDescription
    @test projectTime(mdfFromMem) == testDate
    @test projectUuid(mdfFromMem) == testUuid

    project = defaultMDFv2Project()
    @test !ismissing(projectUuid(project))
    @test !ismissing(projectTime(project))
    @test ismissing(projectName(project))
    @test ismissing(projectDescription(project))

    inMemory = defaultMDFv2InMemory()
    @test !ismissing(projectUuid(inMemory))
    @test !ismissing(projectTime(inMemory))
    @test ismissing(projectName(inMemory))
    @test ismissing(projectDescription(inMemory))
  end
end