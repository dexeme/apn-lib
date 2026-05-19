using Test
using Nemo
using APNLib

include("table1_data.jl")
include("fixtures/table1_luts.jl")

@testset "Kaleyski Table 1 generated LUT fixture" begin
    for case in KALEYSKI_TABLE1_CASES
        key = (case.n, case.id)
        @test haskey(KALEYSKI_TABLE1_GENERATED_LUTS, key)
        @test length(KALEYSKI_TABLE1_GENERATED_LUTS[key]) == 2^case.n
    end
end

@testset "Kaleyski Table 1 implemented permutation counts" begin
    for case in KALEYSKI_TABLE1_CASES
        @testset "n = $(case.n), ID $(case.id)" begin
            lut = KALEYSKI_TABLE1_GENERATED_LUTS[(case.n, case.id)]
            results = reconstruct_external_linear_maps(lut, lut, case.n)

            obtained = length(results)
            expected = kaleyski_table1_expected_permutations(case.n, case.id)

            @info "Kaleyski Table 1 regression case" n = case.n id = case.id obtained expected
            @test obtained == expected
        end
    end
end
