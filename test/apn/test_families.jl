using Test
using APNLib

@testset "APN family recognition" begin
    c4 = APNFunction(monomial_expr(3), Tr(7, monomial_expr(9)))

    c3 = APNFunction(
        8,
        monomial_expr(3),
        monomial_expr(17),
        monomial_expr(48, 18),
        monomial_expr(3, 33),
        monomial_expr(1, 34),
        monomial_expr(48),
    )

    c5 = APNFunction(monomial_expr(3), Tr(39, monomial_expr(9), monomial_expr(18)))

    @test belongs_to_family(c3, :C3)
    @test belongs_to_family(c4, :C4)
    @test belongs_to_family(c5, :C5)

    c3_match = only(filter(match -> match.family == :C3, classify_family(c3)))

    @test c3_match.family == :C3
    @test c3_match.parameters[:m] == 4
    @test c3_match.parameters[:q] == 16
    @test c3_match.parameters[:i] == 1
    @test c3_match.exact == false
end