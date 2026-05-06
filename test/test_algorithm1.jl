using Test
using APNLib

@testset "Algorithm 1 APNsearch Entry Point" begin
    identity1 = reshape([1], 1, 1)

    solutions = APNSearch(
        1,
        identity1,
        identity1,
        max_solutions = 1,
        on_solution = sbox -> nothing,
    )

    @test solutions == [[0, 1]]
end

@testset "Algorithm 1 Function Alias" begin
    identity1 = reshape([1], 1, 1)

    @test APNsearch(1, identity1, identity1, max_solutions = 1, on_solution = sbox -> nothing) == [[0, 1]]
end
