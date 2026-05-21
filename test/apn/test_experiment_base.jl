@testset "Experiment base" begin
    catalogue = Catalogue(
        APNFunction(3, "gold", monomial_expr(3)),
        APNFunction(4, "gold", monomial_expr(3)),
    )

    context = ExperimentContext(
        :example,
        catalogue;
        description = "Example reproducible context",
        source = "unit test",
        root_dir = @__DIR__,
        fixtures = Dict(:values => Dict((3, "gold") => 8, (4, "gold") => 16)),
        expected = Dict(:size => Dict((3, "gold") => 8, (4, "gold") => 16)),
        parameters = Dict(:label => "gold"),
    )

    spec = ExperimentSpec(
        :size,
        (ctx, function_) -> (
            key = catalogue_key(function_),
            obtained = fixture(ctx, :values)[catalogue_key(function_)],
            expected = expected_value(ctx, :size, function_),
            label = parameter(ctx, :label),
        );
        description = "Reads local fixture values",
        columns = [:key, :obtained, :expected, :label],
    )

    result = run_experiment(context, spec, metadata = RunMetadata(env_keys = String[]))

    @test spec.id == :size
    @test spec.columns == [:key, :obtained, :expected, :label]
    @test result.context_id == :example
    @test length(result.rows) == 2
    @test result.rows[1].key == (3, "gold")
    @test result.rows[1].obtained == result.rows[1].expected
    @test result.rows[2].key == (4, "gold")
    @test result.metadata.threads >= 1
end
