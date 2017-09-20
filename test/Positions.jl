@testset "Testing Positions submodule" begin
  shp = [3,3,3]
  fov = [3.0,3.0,3.0]u"mm"
  ctr = [0,0,0]u"mm"
  caG = CartesianGridPositions(shp,fov,ctr)
  @test shape(caG) == shp
  @test fieldOfView(caG) == fov
  @test fieldOfViewCenter(caG) == ctr
  #TODO following tests fail. dont know why
  #@test_throws BoundsError rG[0]
  #@test_throws BoundsError rG[28]
  @test caG[1] == [-1,-1,-1]u"mm"
  @test caG[2] == [0,-1,-1]u"mm"
  @test caG[3] == [1,-1,-1]u"mm"
  @test caG[4] == [-1,0,-1]u"mm"
  @test caG[27] == [1,1,1]u"mm"
  h5open("Positions.h5", "w") do file
    write(file, caG)
  end
  h5open("Positions.h5", "r") do file
    caG1 = CartesianGridPositions(file)
    @test shape(caG1) == shp
    @test fieldOfView(caG1) == fov
    @test fieldOfViewCenter(caG1) == ctr
  end


  chG = ChebyshevGridPositions(shp,fov,ctr)
  @test shape(chG) == shp
  @test fieldOfView(chG) == fov
  @test fieldOfViewCenter(chG) == ctr
  #TODO following tests fail. dont know why
  #@test_throws BoundsError cG[0]
  #@test_throws BoundsError cG[28]
  @test chG[1] ≈ cos(π/6)*3/2*caG[1]
  @test chG[2] ≈ cos(π/6)*3/2*caG[2]
  @test chG[3] ≈ cos(π/6)*3/2*caG[3]
  @test chG[4] ≈ cos(π/6)*3/2*caG[4]
  @test chG[27] ≈ cos(π/6)*3/2*caG[27]
  h5open("Positions.h5", "w") do file
    write(file, chG)
  end
  h5open("Positions.h5", "r") do file
    chG1 = ChebyshevGridPositions(file)
    @test shape(chG1) == shp
    @test fieldOfView(chG1) == fov
    @test fieldOfViewCenter(chG1) == ctr
  end

  mG = MeanderingGridPositions(caG)
  @test shape(mG) == shp
  @test fieldOfView(mG) == fov
  @test fieldOfViewCenter(mG) == ctr
  @test mG[1] == caG[1]
  @test mG[2] == caG[2]
  @test mG[3] == caG[3]
  @test mG[4] == caG[6]
  @test mG[7] == caG[7]
  @test mG[9] == caG[9]
  @test mG[10] == caG[18]
  @test mG[18] == caG[10]
  @test mG[19] == caG[19]
  @test mG[27] == caG[27]
  @test getPermutation(mG) == [1, 2, 3, 6, 5, 4, 7, 8, 9, 18, 17, 16, 13, 14, 15, 12, 11, 10, 19, 20, 21, 24, 23, 22, 25, 26, 27]
  h5open("Positions.h5", "w") do file
    write(file, mG)
  end
  h5open("Positions.h5", "r") do file
    mG1 = MeanderingGridPositions(file)
    @test typeof(mG1.grid) <: CartesianGridPositions
    @test shape(mG1) == shp
    @test fieldOfView(mG1) == fov
    @test fieldOfViewCenter(mG1) == ctr
  end

  positions = [1 2 3 4; 0 1 2 3;-4 -3 -2 -1]u"mm"
  aG1 = ArbitraryPositions(positions)
  @test aG1[1] == [1,0,-4]*u"mm"
  @test aG1[2] == [2,1,-3]u"mm"
  @test aG1[3] == [3,2,-2]u"mm"
  @test aG1[4] == [4,3,-1]u"mm"
  aG2 = ArbitraryPositions(caG)
  @test aG2[1] == caG[1]
  @test aG2[2] == caG[2]
  @test aG2[27] == caG[27]
  h5open("Positions.h5", "w") do file
    write(file, aG2)
  end
  h5open("Positions.h5", "r") do file
    aG3 = ArbitraryPositions(file)
    @test typeof(aG3) <: ArbitraryPositions
    @test aG3.positions == aG2.positions
  end
    
  # the same seed yields the same sequence of points
  seed = UInt32(42)
  N = UInt(3)
  domain = AxisAlignedBox(fov,ctr)
  @test domain.fov == fov
  @test domain.center == ctr
  rG1 = UniformRandomPositions(N,seed,domain)
  @test rG1[1] == [0.09954904813158394,-0.13791259323857274,-1.446939519855107]u"mm"
  @test rG1[2] == [-0.9812009131891462,1.3767776289892044,1.4206979394110573]u"mm"
  @test rG1[3] == [-0.5883911667396526,-0.9692742011014337,1.3707474722677764]u"mm"
  h5open("Positions.h5", "w") do file
    write(file, rG1)
  end
  h5open("Positions.h5", "r") do file
    rG2 = UniformRandomPositions(file)
    @test typeof(rG2) <: UniformRandomPositions{AxisAlignedBox}
    @test rG2.N == N
    @test rG2.seed == seed
    @test rG2.domain.fov == fov
    @test rG2.domain.center == ctr
  end
  radius = 10u"mm"
  domain = Ball(radius,ctr)
  @test domain.radius == radius
  @test domain.center == ctr
  rG3 = UniformRandomPositions(N,seed,domain)
  @test rG3[1] == [-6.715713750009747,0.4103832286623821,-4.525933276650638]u"mm"
  h5open("Positions.h5", "w") do file
    write(file, rG3)
  end
  h5open("Positions.h5", "r") do file
    rG4 = UniformRandomPositions(file)
    @test typeof(rG4) <: UniformRandomPositions{Ball}
    @test rG4.N == N
    @test rG4.seed == seed
    @test rG4.domain.radius == radius
    @test rG4.domain.center == ctr
  end
  #TODO conversion methods dont work. Why?
  #rG = UniformRandomPositions(15,fov,ctr)

  # the following 2 tests fail but should work
  #@test_throws DomainError loadTDesign(pathtosrc,8,1)
  #@test_throws DomainError loadTDesign(pathtosrc,10,1)
  t = 1
  N = 2
  tDesign = loadTDesign(t,N, 5u"mm")
  @test length(tDesign) == N
  @test tDesign.T == t
  @test tDesign.radius == 5u"mm"
  @test any(tDesign.positions .== [1 -1; 0 0; 0 0])
  @test tDesign[1] == [5,0,0]u"mm"
  @test tDesign[2] == [-5,0,0]u"mm"

  @test length(caG) == prod(shp)
  @test length(chG) == prod(shp)
  @test length(mG) == prod(shp)
  @test length(aG1) == size(positions,2)
  @test length(aG2) == prod(shp)

  for (i,p) in enumerate(caG)
    @test p == caG[i]
  end
end
