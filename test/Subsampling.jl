@testset "Positions Subsampling" begin
  @testset "SubsampledPositions" begin
    grid = RegularGridPositions((4, 3), (3.0, 2.0), (0.0, 0.0))  # 12 points
    @test length(grid) == 12

    # Deterministic sampling with a fixed seed
    seed = 0x42
    sub1 = SubsampledPositions(grid, 5; seed = seed)
    sub2 = SubsampledPositions(grid, 5; seed = seed)
    @test length(sub1) == 5
    @test parent(sub1) === grid
    @test parentindices(sub1) == sub1.indices
    @test sub1.indices == sub2.indices  # reproducible
    @test all(1 <= i <= length(grid) for i in sub1.indices)
    @test length(unique(sub1.indices)) == 5
    @test length(collect(sub1)) == length(sub1)

    # getindex should forward to parent at selected indices
    for i in 1:length(sub1)
      @test sub1[i] == grid[sub1.indices[i]]
    end

    # Query parent grid for "metadata"
    @test fieldOfView(sub1) == fieldOfView(grid)
    @test fieldOfViewCenter(sub1) == fieldOfViewCenter(grid)
    @test shape(sub1) == shape(grid)
    @test spacing(sub1) == spacing(grid)

    # User-specified indices, ranges
    sub3 = SubsampledPositions(grid, 1:2:4)
    @test parentindices(sub3) == collect(1:2:4)
    @test length(sub3) == 2
    @test sub3[2] == grid[(1:2:4)[end]]
    # User-specified indices, preserve order
    sub4 = SubsampledPositions(grid, reverse(1:length(grid)))
    for i in 1:length(sub4)
      @test sub4[i] == grid[12 - (i - 1)]
    end

    # View-constructor
    sub5 = view(grid, 1:2:4)
    @test parentindices(sub5) == parentindices(sub3)

    # Positions constructor
    sub6 = SubsampledPositions(grid, hcat(collect(sub1)...))
    @test parentindices(sub6) == parentindices(sub1)
    sub6a = SubsampledPositions(grid.shape, grid.fov, grid.center, hcat(collect(sub1)...))
    @test parentindices(sub6a) == parentindices(sub1)

    # Different seeds should typically produce different selections
    sub_diff = SubsampledPositions(grid, 5; seed = 0xdeadbeef)
    @test sub_diff.indices != sub1.indices

    # Factor-based constructor: just check the resulting length
    subf = SubsampledPositions(grid, 0.25; seed = seed)
    @test length(subf) == round(Int, length(grid) * 0.25)

    @testset "Edge cases" begin
      grid = RegularGridPositions((1, 1), (0.0, 0.0), (0.0, 0.0))  # single point
      @test length(grid) == 1
      sub = SubsampledPositions(grid, 1; seed = 0x1)
      @test sub.indices == [1]
      @test sub[1] == grid[1]

      # Empty subsample
      sub_empty = SubsampledPositions(grid, 0; seed = 0x1)
      @test length(sub_empty) == 0
      @test isempty(parentindices(sub_empty))

      # Invalid arguments, too many indices
      @test_throws ArgumentError SubsampledPositions(grid, length(grid) + 1)
      # Invalid arguments, invalid indices
      @test_throws ArgumentError SubsampledPositions(grid, [length(grid) + 1])
      # Invalid arguments, incorrect position matrix size, > D
      @test_throws ArgumentError SubsampledPositions(grid, fill(0.0, length(first(grid)) + 1, 1))
      # Invalid arguments, incorrect position matrix size, < D 
      @test_throws ArgumentError SubsampledPositions(grid, fill(0.0, length(first(grid)) - 1, 1))
      # Invalid arguments, no matching position found
      @test_throws ArgumentError SubsampledPositions(grid, fill(1.0, length(first(grid)), 1))
    end

    @testset "Nesting" begin
      # Test "SM" with bg measurement:
      grid = RegularGridPositions((3, 3, 3), (3.0, 2.0, 3.0), (0.0, 0.0, 0.0)) # 27 points
      csgrid = SubsampledPositions(grid, 1:2:length(grid))
      # BG Meas at the start, middle and end
      bgInd = Int64.(collect(range(1, length(csgrid) + 3, length = 3)))
      bgPos = [10.0, 10.0, 10.0]
      bggrid = BreakpointPositions(csgrid, bgInd, bgPos)
      @test length(bggrid) == length(csgrid) + length(bgInd)
      @test shape(bggrid) == shape(grid)
      @test fieldOfView(bggrid) == fieldOfView(grid)
      @test fieldOfViewCenter(bggrid) == fieldOfViewCenter(grid)
      @test bggrid[1] == bgPos
      @test bggrid[2 + div(length(csgrid), 2)] == bgPos
      @test bggrid[length(bggrid)] == bgPos

      # Test "SM" with bg measurements and sorting
      grid = RegularGridPositions((3, 3, 3), (3.0, 2.0, 3.0), (0.0, 0.0, 0.0)) # 27 points
      csgrid = SubsampledPositions(grid, 1:2:length(grid))
      sortedgrid = SortedPositions(csgrid)
      # BG Meas at the start, middle and end
      bgInd = Int64.(collect(range(1, length(csgrid) + 3, length = 3)))
      bgPos = [10.0, 10.0, 10.0]
      bggrid = BreakpointPositions(csgrid, bgInd, bgPos)
      @test length(bggrid) == length(csgrid) + length(bgInd)
      @test shape(bggrid) == shape(grid)
      @test fieldOfView(bggrid) == fieldOfView(grid)
      @test fieldOfViewCenter(bggrid) == fieldOfViewCenter(grid)
      @test bggrid[1] == bgPos
      @test bggrid[2 + div(length(csgrid), 2)] == bgPos
      @test bggrid[length(bggrid)] == bgPos
    end
  end

end