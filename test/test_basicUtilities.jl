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
