using Test
using Nemo
using APNLib

include("KaleyskiExperiments.jl")
using .KaleyskiExperiments

@testset "Kaleyski Table 1 generated LUT fixture" begin
    context = KALEYSKI_TABLE1_CONTEXT

    @test context.catalogue === KaleyskiExperiments.KALEYSKI_TABLE1_CATALOGUE
    @test haskey(context.fixtures, :luts)
    @test haskey(context.expected, :external_linear_permutation_count)

    for case in KaleyskiExperiments.KALEYSKI_TABLE1_CASES
        key = (case.n, case.id)
        @test haskey(fixture(context, :luts), key)
        @test length(fixture(context, :luts)[key]) == 2^case.n
    end
end

@testset "Kaleyski Table 1 implemented permutation counts" begin
    for case in KaleyskiExperiments.KALEYSKI_TABLE1_CASES
        @testset "n = $(case.n), ID $(case.id)" begin
            context = kaleyski_table1_selected_context([case])
            result = run_kaleyski_table1_experiment(context = context)
            row = only(result.rows)

            @info "Kaleyski Table 1 regression case" n = row.n id = row.id obtained = row.found expected = row.expected
            @test row.ok
        end
    end
end
