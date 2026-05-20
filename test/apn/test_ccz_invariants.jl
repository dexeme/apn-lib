using Test
using Nemo
using APNLib

@testset "CCZ invariants" begin
    identity_lut = [0, 1, 2, 3]

    @test gamma_rank(identity_lut, 2) isa Int
    @test delta_rank(identity_lut, 2) isa Int

    n = 5
    gold_polynomial = APNFunction(monomial_expr(3))
    gold_lut = univariate_to_lut(gold_polynomial, n)

    @test is_apn(gold_lut)

    gamma = gamma_rank(gold_lut, n)
    delta = delta_rank(gold_lut, n)

    @test gamma == 330
    @test delta == 42
    @test gamma_rank(gold_polynomial, n) == gamma
    @test delta_rank(gold_polynomial, n) == delta

    dimensioned_gold_polynomial = APNFunction(n, monomial_expr(3))

    @test gamma_rank(dimensioned_gold_polynomial, n) == gamma
    @test delta_rank(dimensioned_gold_polynomial, n) == delta
end
