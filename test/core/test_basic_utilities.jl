using Test
using Nemo
using APNLib

@testset "Basic Utility Functions" begin
    F = gf2()

    identity3 = identity_matrix(F, 3)

    @test check_square(identity3) == true
    @test check_same_size(identity3, identity3) == true
    @test check_same_field(identity3, identity3) == true
    @test check_compatible_pair(identity3, identity3) == true

    @test int_to_bits(5, 3) == [1, 0, 1]
    @test column_vector_to_int(int_to_column_vector(F, 5, 3)) == 5
    @test length(collect(each_column_vector(F, 3))) == 8
end

@testset "Core Validation Helpers" begin
    F = gf2()
    K = GF(2, 2, "a")
    identity2 = identity_matrix(F, 2)
    identity3 = identity_matrix(F, 3)
    extension_identity2 = identity_matrix(K, 2)

    @test space_size(0) == 1
    @test space_size(4) == 16
    @test_throws ErrorException space_size(-1)

    @test check_length([1, 2, 3], 3) == true
    @test_throws ErrorException check_length([1, 2], 3)
    @test check_space_length(zeros(Int, 8), 3) == true
    @test_throws ErrorException check_space_length(zeros(Int, 7), 3)

    @test check_integer_range(2, 1, 3) == true
    @test_throws ErrorException check_integer_range(4, 1, 3)
    @test check_space_value(7, 3) == true
    @test_throws ErrorException check_space_value(8, 3)

    @test check_n_by_n_matrix(identity3, 3) == true
    @test_throws ErrorException check_n_by_n_matrix(identity3, 2)
    @test check_gf2_matrix(identity2, 2) == true
    @test_throws ErrorException check_gf2_matrix(extension_identity2)

    @test check_same_size(identity2, identity2) == true
    @test_throws ErrorException check_same_size(identity2, identity3)
    @test_throws ErrorException check_same_field(identity2, extension_identity2)
    @test_throws ErrorException check_compatible_pair(identity2, identity3)

    @test check_lut_values([0, 1, 2, 3], 2) == true
    @test_throws ErrorException check_lut_values([0, 1, 2], 2)
    @test_throws ErrorException check_lut_values([0, 1, 2, 4], 2)
    @test check_sbox_ddt_sizes([0, 1, 2, 3], zeros(Int, 4, 4)) == true
    @test_throws ErrorException check_sbox_ddt_sizes([0, 1, 2, 3], zeros(Int, 4, 3))
end

@testset "Finite Field Bit Conversions" begin
    F = gf2()

    @test int_to_bits(0, 4) == [0, 0, 0, 0]
    @test int_to_bits(6, 4) == [0, 1, 1, 0]
    @test int_to_column_vector(F, 6, 4) == matrix(F, 4, 1, [0, 1, 1, 0])
    @test column_vector_to_int(matrix(F, 4, 1, [0, 1, 1, 0])) == 6
    @test [column_vector_to_int(v) for v in each_column_vector(F, 2)] == collect(0:3)
end

@testset "Core Linear Algebra Helpers" begin
    F = gf2()
    identity2 = identity_matrix(F, 2)
    swap = matrix(F, 2, 2, [0, 1, 1, 0])
    block = APNLib.block_diagonal_gf2([identity2, swap])

    @test nrows(block) == 4
    @test ncols(block) == 4
    @test block[1:2, 1:2] == identity2
    @test block[3:4, 3:4] == swap
    @test matrix_multiplicative_order(identity2) == 1
    @test matrix_multiplicative_order(swap) == 2
    @test_throws ErrorException APNLib.block_diagonal_gf2(FqMatrix[])
end
