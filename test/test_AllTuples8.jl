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

    @testset "Proposition 5" begin
        # 3, 7, 11, 12, 17, 21 and 24 should be filtered by Proposition 5
        classes_filtered_by_prop_5 = [3, 7, 11, 12, 17, 21, 24]
        for i in classes_filtered_by_prop_5
            A, B = precomputed_tuple_matrices(8, i)
            @test filtro_proposicao_5(A, B, 8) == false
        end

        # 1, 2, 4:6, 8:10, 13:16, 18:20, 22, 23, 25:32 should pass Proposition 5
        classes_that_pass_prop_5 = [1, 2, 4, 5, 6, 8, 9, 10, 13, 14, 15, 16, 18, 19, 20, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32]
        for i in classes_that_pass_prop_5
            A, B = precomputed_tuple_matrices(8, i)
            @test filtro_proposicao_5(A, B, 8) == true
        end
    end
end
