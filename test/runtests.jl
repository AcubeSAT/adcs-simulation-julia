using ADCSSims
using Test
using LinearAlgebra

@testset "Frame alignment" begin
    S, N = normalize(rand(3)), normalize(rand(3))
    Y, Z = [0,1,0], [0,0,-1]
    q = align_frame_with_vector(S, N, Z, Y)
    Y_new = rotvec(Y,q)
    Z_new = rotvec(Z,q)
    @test Z_new ≈ S && dot(Y_new, Y) ≥ 0.0
end