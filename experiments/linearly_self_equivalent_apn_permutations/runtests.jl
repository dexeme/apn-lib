using Test
using Nemo
using APNLib

include("tuple_generation.jl")
include("test_proposition4.jl")
include("test_proposition5.jl")

@testset "Linearly self-equivalent APN permutation tuple experiment" begin
    F = gf2()
    identity3 = identity_matrix(F, 3)

    @test check_order_space(identity3, 2) == 3.0
    @test is_permutation_tuple(identity3, identity3) == true
    @test length(gen_permutation_tuples(6)) == 17

    # T7 = gen_permutation_tuples(7)
    # @test length(T7) == 27

    # T8 = gen_permutation_tuples(8)
    # @test length(T8) == 32
end
