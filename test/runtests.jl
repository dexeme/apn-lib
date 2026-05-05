using Test
using Nemo
using APNLib

@testset "Basic utilities" begin
    F = gf2()

    A = identity_matrix(F, 3)
    B = identity_matrix(F, 3)

    @test check_order_space(A, 2) == 3.0
    @test matrix_to_sbox(A) == collect(0:7)
    @test is_permutation_tuple(A, B) == true
end

@testset "RCF blocks" begin
    blocks = blocks_for_rcf(3)

    @test length(blocks) == 4
end

@testset "RCFs" begin
    @test length(get_rcfs(2)) == 3
    @test length(get_rcfs(3)) == 6
    @test length(get_rcfs(4)) == 14
end

@testset "Generation counts" begin
    println("a")
    T6 = gen_permutation_tuples(6)
    @test length(T6) == 17

    # T7 = gen_permutation_tuples(7)
    # @test length(T7) == 27
    
    # T8 = gen_permutation_tuples(8)
    # @test length(T8) == 32
end