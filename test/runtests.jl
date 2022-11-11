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

    err = []
    try
        El0 = Ellipsoid(-1*P0, c0)
    catch e
        err = e
    end

    @test err isa Exception
    @test sprint(showerror, err) == "P must be a positive definite matrix"


end


@testset "EllipsoidIncluded" begin
    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]

    c0 = [1.6; 1.4]
    P0 = [0.4 -0.1;
         -0.1 0.5]

    El = Ellipsoid(P, c)
    El0 = Ellipsoid(P0, c0)

    @test El ∈ El0
    @test 6*El ∈ El0
    @test El0 ∈ El0
    @test 0.5*El0 ∈ El0
    @test 0.5*El ∈ El0
    
    ell_ast, β_ast =  EllipsoidInclusion.get_ℓ_ast(El, El0)
    @test all(isapprox.(EllipsoidInclusion.get_ℓ_ast(El, El0*(-ell_ast)), (-1.0, -β_ast/ell_ast); atol=1e-4))

    @test El ∈ El0*(-ell_ast)*1.00001
    @test El ∉ El0*(-ell_ast)*0.99999
end


@testset "EllipsoidNotIncluded" begin
    c = [1.6; 1.4]
    P = [0.4 -0.1;
         -0.1 0.5]
    
    c0 = [1.5; 1.5]
    P0 = [4.0 0.5;       
         0.5 6.0]

    El = Ellipsoid(P, c)
    El0 = Ellipsoid(P0, c0)

    @test El ∉ El0
    @test El*0.05 ∉ El0
    @test El*10 ∉ El0
    @test El0*1.01 ∉ El0
    @test El ∉ El0*0.9
    
    ell_ast, β_ast =  EllipsoidInclusion.get_ℓ_ast(El, El0)
    @test all(isapprox.(EllipsoidInclusion.get_ℓ_ast(El, El0*(-ell_ast)), (-1.0, -β_ast/ell_ast); atol=1e-4))

    @test El ∈ El0*(-ell_ast)*1.00001
    @test El ∉ El0*(-ell_ast)*0.99999
end



end # module
