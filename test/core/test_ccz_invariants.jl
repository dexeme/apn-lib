using Test
using Nemo
using APNLib

@testset "CCZ invariants" begin
    identity_lut = [0, 1, 2, 3]

    @test gamma_rank(identity_lut, 2) isa Int
    @test delta_rank(identity_lut, 2) isa Int

    n = 5
    field = GF(2, n, "g")
    R, x = Nemo.polynomial_ring(field, "x")
    gold_lut = univariate_to_lut(x^3, n)

    @test is_apn(gold_lut)

    gamma = gamma_rank(gold_lut, n)
    delta = delta_rank(gold_lut, n)

    @test gamma == 330
    @test delta == 42
end