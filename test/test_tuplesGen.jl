using Test
using Nemo
using APNLib

@testset "Tuple Generation And Precomputed Constants" begin
    F = gf2()
    identity3 = identity_matrix(F, 3)

    @test check_order_space(identity3, 2) == 3.0
    @test matrix_to_sbox(identity3) == collect(0:7)
    @test permutation_cycle_structure([0, 2, 1, 3]) == [1, 1, 2]
    @test matrix_cycle_structure(identity3) == ones(Int, 8)
    @test same_cycle_structure(identity3, identity3) == true
    @test is_permutation_tuple(identity3, identity3) == true

    rows = load_precomputed_tuple_constants(6)
    lut_A, lut_B = precomputed_tuple_sboxes(6, 1)

    @test length(rows) == 17
    @test precomputed_tuple_row(6, 17) == rows[17]
    @test length(precomputed_tuple_row(6, 1)) == 128
    @test length(lut_A) == 64
    @test length(lut_B) == 64
    @test vcat(lut_B, lut_A) == rows[1]
end

@testset "Tuple Matrix Constant Extraction" begin
    row = precomputed_tuple_row(6, 1)
    A_from_lut, B_from_lut = extract_matrices(row, 6)
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
end

@testset "Permutation Tuple Generation Counts" begin
    tuples6 = gen_permutation_tuples(6)
    @test length(tuples6) == 17

    # T7 = gen_permutation_tuples(7)
    # @test length(T7) == 27

    # T8 = gen_permutation_tuples(8)
    # @test length(T8) == 32
end
