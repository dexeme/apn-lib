using Test
using Nemo
using APNLib

@testset "Gold function representation conversions" begin
    n = 3
    gold_polynomial = APNFunction(monomial_expr(3))

    gold_graph = univariate_to_graph(gold_polynomial, n)
    gold_lut = graph_to_lut(gold_graph, n)
    gold_anf = graph_to_anf(gold_graph, n)
    recovered_polynomial, _ = anf_to_univariate(gold_anf)

    @test length(gold_graph) == 2^n
    @test anf_to_lut(gold_anf) == gold_lut
    @test univariate_to_lut(recovered_polynomial, n) == gold_lut
    @test univariate_to_lut(gold_polynomial, n) == gold_lut
end

@testset "Graph and LUT validation" begin
    n = 2

    @test lut_to_graph([0, 1, 2, 3], n) == [(0, 0), (1, 1), (2, 2), (3, 3)]
    @test graph_to_lut([(2, 3), (0, 1), (3, 0), (1, 2)], n) == [1, 2, 3, 0]
    @test_throws ErrorException graph_to_lut([(0, 0), (0, 1), (2, 2), (3, 3)], n)
    @test_throws ErrorException graph_to_lut([(0, 0), (1, 1), (2, 2)], n)
    @test_throws ErrorException lut_to_anf([0, 1, 2, 4], n)
end
