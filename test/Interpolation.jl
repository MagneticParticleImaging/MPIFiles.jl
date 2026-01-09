@testset "Testing Interpolation submodule" begin
  
  grid = RegularGridPositions([5,5,10],[1.0,1.0,1.0],[0.0,0.0,0.0])
  A = randn(5,5,10)
  @test all(A .≈ MPIFiles.interpolate(A, grid, grid))

  origin = RegularGridPositions([3,3,1],[3.0,3.0,1.0],[0.0,0.0,0.0])
  target = RegularGridPositions([3,3,1],[3.0,3.0,1.0],[0.0,0.0,0.0], [-1,-1,1])
  A = reshape(1.0:9.0,Tuple(origin.shape))
  @test all(reverse(A) .≈ MPIFiles.interpolate(A, origin, target))

  A = Float64.(reshape([x+y-2 for y in 1:3 for x in 1:3],Tuple(origin.shape)))
  origin = RegularGridPositions([3,3,1],[3.0,3.0,1.0],[0.0,0.0,0.0])
  target = axesToRegularGridPositions(-1:0.5:1, -1:0.5:1, 0.:0.)
  @test all(A .≈ MPIFiles.interpolate(A, origin, target)[1:2:end, 1:2:end, 1:end])
end