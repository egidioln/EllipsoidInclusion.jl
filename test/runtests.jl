module TestMain
using Test
using EllipsoidInclusion

sleep(0.1) # used for good printing
println("Started test")

@testset "EllipsoidCreation" begin
    n = 10
    aux = randn(n,n)
    P0 = aux'aux
    c0 = randn(n)
    
    
    El0 = Ellipsoid(P0, c0)
    
    @test EllipsoidInclusion.get_center(El0) == c0
    @test c0 ∈ El0

    end

end # module
