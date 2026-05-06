# test_AllTuples6 Tests

using Test
using Nemo
using APNLib

@testset "test_AllTuples6 Tests" begin
    @testset "Proposition 4" begin
        # 6, 9, 13, 16 and 17 should be filtered by Proposition 4
        classes_filtered_by_prop_4 = [6, 9, 13, 16, 17]
        for i in classes_filtered_by_prop_4
            A, B = precomputed_tuple_matrices(6, i)
            @test filtro_proposicao_4(A, B, 6) == false
        end

        # 1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 14 and 15 should pass Proposition 4
        classes_that_pass_prop_4 = [1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 14, 15]
        for i in classes_that_pass_prop_4
            A, B = precomputed_tuple_matrices(6, i)
            @test filtro_proposicao_4(A, B, 6) == true
        end




    end
end
