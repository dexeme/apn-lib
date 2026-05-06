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
    @test filtro_proposicao_4(A, B, 3) == true
    identity3_array = [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test filtro_proposicao_4(identity3_array, identity3_array, 3) == true
end

@testset "Proposition 4 filter" begin
    F = gf2()

    identity3 = identity_matrix(F, 3)
    transvection = matrix(F, 3, 3, [
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
    ])
    cycle = matrix(F, 3, 3, [
        0, 1, 0,
        0, 0, 1,
        1, 0, 0,
    ])

    @test filtro_proposicao_4(transvection, transvection, 3) == false
    @test filtro_proposicao_4(identity3, cycle, 3) == false
end

@testset "Precomputed tuple constants" begin
    rows = load_precomputed_tuple_constants(6)
    lut_A, lut_B = precomputed_tuple_sboxes(6, 1)

    @test length(rows) == 17
    @test precomputed_tuple_row(6, 17) == rows[17]
    @test length(precomputed_tuple_row(6, 1)) == 128
    @test length(lut_A) == 64
    @test length(lut_B) == 64
    @test vcat(lut_B, lut_A) == rows[1]
end

@testset "Precomputed tuple matrix constants" begin
    row = precomputed_tuple_row(6, 1)
    A_from_lut, B_from_lut = extrair_matrizes(row, 6)
    A_from_file, B_from_file = precomputed_tuple_matrices(6, 1)

    expected_A = [
        0 0 0 0 0 1
        1 0 0 0 0 0
        0 1 0 0 0 1
        0 0 1 0 0 1
        0 0 0 1 0 0
        0 0 0 0 1 0
    ]
    expected_B = [
        0 0 0 0 0 1
        1 0 0 0 0 1
        0 1 0 0 0 0
        0 0 1 0 0 1
        0 0 0 1 0 1
        0 0 0 0 1 1
    ]

    @test A_from_lut == expected_A
    @test B_from_lut == expected_B
    @test A_from_file == A_from_lut
    @test B_from_file == B_from_lut
    @test filtro_proposicao_4(A_from_file, B_from_file, 6) == true
end

@testset "Dynamic APN search utilities" begin
    sbox = fill(-1, 4)
    ddt = zeros(Int, 4, 4)

    sbox[0 + 1] = 0
    @test addDDTInformation(0, sbox, ddt) == true
    @test all(ddt .== 0)

    sbox[3 + 1] = 3
    @test addDDTInformation(3, sbox, ddt) == true
    @test ddt[3 + 1, 3 + 1] == 2
    @test removeDDTInformation(3, sbox, ddt) == true
    @test all(ddt .== 0)

    @test isComplete(sbox) == false
    @test nextFreePosition(sbox) == 1

    A = reshape([1], 1, 1)
    solutions = APNSearch(1, A, A, max_solutions = 1, on_solution = s -> nothing)
    @test solutions == [[0, 1]]
end

include("test_AllTuples6.jl")
include("test_AllTuples7.jl")
include("test_AllTuples8.jl")

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
