# test_AllTuples8 Tests

using Test
using Nemo
using APNLib

@testset "test_AllTuples8 Tests" begin
    @testset "Proposition 4" begin
        # 13, 18, 19, 20, 21, 25, 26, 28, 32 should be filtered by Proposition 4
        classes_filtered_by_prop_4 = [13, 18, 19, 20, 21, 25, 26, 28, 32]
        for i in classes_filtered_by_prop_4
            A, B = precomputed_tuple_matrices(8, i)
            @test filtro_proposicao_4(A, B, 8) == false
        end

        # 1:12, 14:17, 22:24, 27, 29:31 should pass Proposition 4
        classes_that_pass_prop_4 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 22, 23, 24, 27, 29, 30, 31]
        for i in classes_that_pass_prop_4
            A, B = precomputed_tuple_matrices(8, i)
            @test filtro_proposicao_4(A, B, 8) == true
        end




    end
end
