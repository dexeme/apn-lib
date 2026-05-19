using Test
using Nemo
using APNLib

@testset "EA External Linear Reconstruction" begin
    n = 6
    field = GF(2, n, "a")
    R, z = Nemo.polynomial_ring(field, "x")
    alpha = gen(field)

    gold_lut = univariate_to_lut(z^3, n)
    gold_partition = partition_by_multiplicity(gold_lut, n)
    @test sort(length.(gold_partition.blocks)) == [1, 21, 42]

    gold_results = reconstruct_external_linear_maps(gold_lut, gold_lut, n, parallel = false)
    @test length(gold_results) == 1008

    gold_parallel_results = reconstruct_external_linear_maps(gold_lut, gold_lut, n, parallel = true)
    @test sort(gold_parallel_results) == sort(gold_results)

    restricted_lut = univariate_to_lut(z^3 + alpha^11 * z^6 + alpha * z^9, n)
    restricted_partition = partition_by_multiplicity(restricted_lut, n)
    @test sort(length.(restricted_partition.blocks)) == [1, 21, 42]

    restricted_results = reconstruct_external_linear_maps(restricted_lut, restricted_lut, n)
    @test length(restricted_results) == 336
end
