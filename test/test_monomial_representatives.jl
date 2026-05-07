using Test
using APNLib

@testset "Monomial representative search" begin
    candidates = N7_MONOMIAL_REPRESENTATIVE_EXPONENTS

    @test compute_expected_monomial_representatives(7, 1, candidate_exponents = candidates) ==
          Dict(1 => [5, 9, 63, 78, 85, 88])

    @test compute_expected_monomial_representatives(7, [1, 4, 2], candidate_exponents = candidates) ==
          Dict(
              1 => [5, 9, 63, 78, 85, 88],
              2 => Int[],
              4 => [63],
          )

    expected_search_classes = [1, 4, 5, 7, 8, 9, 10, 16, 17, 18, 20, 22, 23, 25, 26]
    @test search_classes_from_saved_file(7) == expected_search_classes
    all_results = compute_expected_monomial_representatives(7, "all", candidate_exponents = candidates)
    @test sort(collect(keys(all_results))) == expected_search_classes
    @test all_results[1] == [5, 9, 63, 78, 85, 88]
    @test all_results[4] == [63]
end
