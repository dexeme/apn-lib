using Test
using Nemo
using APNLib

@testset "Gold function representation conversions" begin
    n = 3
    field = GF(2, n, "g")
    polynomial_ring, x = Nemo.polynomial_ring(field, "x")
    gold_polynomial = x^3

    gold_graph = univariate_to_graph(gold_polynomial, n)
    gold_lut = graph_to_lut(gold_graph, n)
    gold_anf = graph_to_anf(gold_graph, n)
    recovered_polynomial, _ = anf_to_univariate(gold_anf)

    @test length(gold_graph) == 2^n
    @test anf_to_lut(gold_anf) == gold_lut
    @test univariate_to_lut(recovered_polynomial, n) == gold_lut
    @test univariate_to_lut(gold_polynomial, n) == gold_lut
end

