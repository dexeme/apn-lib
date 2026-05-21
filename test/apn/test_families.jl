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

@testset "APN equivalence catalog" begin
    n = 3
    gold = APNFunction(n, monomial_expr(3))
    gold_lut = apn_to_lut(gold)

    representatives = family_representatives(n)
    gold_representative = only(filter(representative -> representative.id == :gold_i1, representatives))

    @test haskey(APN_REPRESENTATIVES_BY_DIMENSION, n)
    @test representative_lut(gold_representative) == gold_lut
    @test any(representative -> representative.id == :c4, representatives)
    @test any(representative -> representative.id == :c5, representatives)
    @test any(representative -> representative.id == :c6, representatives)

    matches = find_equivalences(n, gold_lut, representatives = [gold_representative])

    @test length(matches) == 1
    @test matches[1].representative.id == :gold_i1
    @test matches[1].equivalence isa EAEquivalence

    @test_throws ArgumentError find_equivalences(n, gold_lut, max_external_maps_per_representative = 0)

    tuple_lut = precomputed_tuple_search(7, 9)
    tuple_matches = find_equivalences(7, tuple_lut)

    @test length(tuple_matches) == 1
    @test tuple_matches[1].representative.id == :tuple_search_9
end
