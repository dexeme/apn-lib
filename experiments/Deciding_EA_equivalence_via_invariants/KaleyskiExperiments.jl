module KaleyskiExperiments

using APNLib

include("table1_data.jl")
include(joinpath("fixtures", "table1_luts.jl"))

export KALEYSKI_TABLE1_CONTEXT,
       KALEYSKI_TABLE1_PERMUTATION_COUNT,
       kaleyski_table1_context,
       kaleyski_table1_permutation_count_spec,
       kaleyski_table1_lut,
       kaleyski_table1_selected_context,
       run_kaleyski_table1_experiment

function kaleyski_table1_context(; catalogue::Catalogue = KALEYSKI_TABLE1_CATALOGUE)
    return ExperimentContext(
        :kaleyski_table1,
        catalogue;
        description = "Reproduction context for Table 1 of Kaleyski, Deciding EA-equivalence via invariants.",
        source = "Kaleyski, Deciding EA-equivalence via invariants; test functions from Edel and Pott's APN catalogue.",
        root_dir = @__DIR__,
        artifacts_dir = joinpath(@__DIR__, "artifacts"),
        fixtures = Dict(:luts => KALEYSKI_TABLE1_GENERATED_LUTS),
        expected = Dict(:external_linear_permutation_count => KALEYSKI_TABLE1_EXPECTED),
        metadata = Dict(
            :case_order => KALEYSKI_TABLE1_CASES,
            :generated_cases => KALEYSKI_TABLE1_GENERATION_CASES,
            :equation_ids => KALEYSKI_TABLE1_EQUATION_IDS,
        ),
    )
end

const KALEYSKI_TABLE1_CONTEXT = kaleyski_table1_context()

function kaleyski_table1_selected_context(cases)
    functions_by_key = Dict(catalogue_key(function_) => function_
                            for function_ in KALEYSKI_TABLE1_FUNCTIONS)
    functions = [functions_by_key[(case.n, case.id)] for case in cases]
    return kaleyski_table1_context(catalogue = Catalogue(functions))
end

function kaleyski_table1_lut(context::ExperimentContext, function_::APNFunction)
    return fixture(context, :luts)[catalogue_key(function_)]
end

function _env_log_level()
    return Symbol(get(ENV, "APNLIB_LOG_LEVEL", "quiet"))
end

function _env_parallel()
    value = lowercase(get(ENV, "APNLIB_PARALLEL", "auto"))
    value == "auto" && return Threads.nthreads() > 1
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("APNLIB_PARALLEL must be auto, true, or false")
end

function _permutation_count_row(context::ExperimentContext, function_::APNFunction)
    parallel = _env_parallel()
    log_level = _env_log_level()
    found = 0

    elapsed = @elapsed begin
        lut = kaleyski_table1_lut(context, function_)
        found = length(reconstruct_external_linear_maps(lut,
                                                        lut,
                                                        function_.n,
                                                        parallel = parallel,
                                                        log_level = log_level))
    end

    expected = expected_value(context, :external_linear_permutation_count, function_)

    return (
        n = function_.n,
        id = function_.id,
        equation = kaleyski_table1_equation_id(function_.n, function_.id),
        time = elapsed,
        found = found,
        expected = expected,
        parallel = parallel,
        log_level = log_level,
        ok = found == expected,
    )
end

function kaleyski_table1_permutation_count_spec()
    return ExperimentSpec(
        :external_linear_permutation_count,
        _permutation_count_row;
        description = "Number of external linear permutations reconstructed for each Table 1 function.",
        columns = [:n, :id, :equation, :time, :found, :expected, :parallel, :log_level, :ok],
    )
end

const KALEYSKI_TABLE1_PERMUTATION_COUNT = kaleyski_table1_permutation_count_spec()

function run_kaleyski_table1_experiment(; context::ExperimentContext = KALEYSKI_TABLE1_CONTEXT,
                                        spec::ExperimentSpec = KALEYSKI_TABLE1_PERMUTATION_COUNT)
    return run_experiment(
        context,
        spec;
        metadata = RunMetadata(env_keys = [
            "APNLIB_PARALLEL",
            "APNLIB_LOG_LEVEL",
            "APNLIB_KALEYSKI_TABLE1_DIMENSIONS",
            "APNLIB_KALEYSKI_TABLE1_IDS",
        ]),
    )
end

end
