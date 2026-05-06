# test_AllTuples7 Tests

using Test
using Nemo
using APNLib

@testset "test_AllTuples7 Tests" begin
    @testset "Proposition 4" begin
        # 12, 13, 14, 15, 20, 21, 25, 27 should be filtered by Proposition 4
        classes_filtered_by_prop_4 = [12, 13, 14, 15, 20, 21, 25, 27]
        for i in classes_filtered_by_prop_4
            A, B = precomputed_tuple_matrices(7, i)
            @test filtro_proposicao_4(A, B, 7) == false
        end

        # 1:11, 16:19, 22:24, 26 should pass Proposition 4
        classes_that_pass_prop_4 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 22, 23, 24, 26]
        for i in classes_that_pass_prop_4
            A, B = precomputed_tuple_matrices(7, i)
            @test filtro_proposicao_4(A, B, 7) == true
        end




    end
end
