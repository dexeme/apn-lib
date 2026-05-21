using Nemo
using APNLib

isdefined(@__MODULE__, :algorithm_1) || include(joinpath(@__DIR__, "algorithm_1.jl"))
isdefined(@__MODULE__, :algorithm_2) || include(joinpath(@__DIR__, "algorithm_2.jl"))

# Data from Table 1 of Kaleyski, "Deciding EA-equivalence via invariants",
# with test functions taken from Edel and Pott's APN catalogue.
#@article{Edel:200900,
# doi = {10.3934/amc.2009.3.59},
# author = {Yves Edel and Alexander Pott},
# title = {{A new almost perfect nonlinear function which is not quadratic}},
# year = 2009,
# journal = {{Advances in Mathematics of Communications}},
# volume = 3,
# number = 1,
# pages = {59--81},
#}
# Keep this file as the readable source of truth for Table 1. The LUTs are
# generated pointwise from these definitions when the experiment module loads.
# In particular, reference(...) means "use the already generated LUT value for
# this input and add the perturbation".

const TABLE_1_EXPECTED = Dict{Tuple{Int, String}, Int}(
    (6, "1.1") => 1008,
    (6, "1.2") => 336,
    (6, "2.1") => 10,
    (6, "2.2") => 336,
    (6, "2.3") => 1008,
    (6, "2.4") => 8,
    (6, "2.5") => 60,
    (6, "2.6") => 8,
    (6, "2.7") => 10,
    (6, "2.8") => 8,
    (6, "2.9") => 8,
    (6, "2.10") => 8,
    (6, "2.11") => 8,
    (6, "2.12") => 48,
    (8, "1.1") => 680,
    (8, "1.2") => 8,    # Kaleyski alocou na linha 1.3
    (8, "1.3") => 4,    # Kaleyski alocou na linha 1.5
    (8, "1.4") => 1,    # A 1ª 'new', Kaleyski alocou na linha 1.7
    (8, "1.5") => 4,    # A 2ª 'new', Kaleyski alocou na linha 1.8
    (8, "1.6") => 4,    # A 3ª 'new', Kaleyski alocou na linha 1.9
    (8, "1.7") => 2,    # A 4ª 'new', Kaleyski alocou na linha 1.10
    (8, "1.8") => 4,    # A 5ª 'new', Kaleyski alocou na linha 1.11
    (8, "1.9") => 4,    # A 6ª 'new', Kaleyski alocou na linha 1.12
    (8, "1.10") => 2,   # A 7ª 'new', Kaleyski alocou na linha 1.13
    (8, "1.11") => 8,   # Kaleyski alocou na linha 1.4
    (8, "1.12") => 4,   # Kaleyski alocou na linha 1.6
    (8, "1.13") => 2,   # A 8ª 'new', Kaleyski alocou na linha 1.14
    (8, "1.14") => 1,   # A 9ª 'new', Kaleyski alocou na linha 1.15
    (8, "1.15") => 680, # Kaleyski alocou na linha 1.2
    (8, "1.16") => 2,   # A 10ª 'new', Kaleyski alocou na linha 1.16
    (8, "1.17") => 2,   # A 11ª 'new', Kaleyski alocou na linha 1.17
    (8, "2.1") => 360,
    (8, "3.1") => 4,
    (8, "4.1") => 16,
    (8, "5.1") => 8,
    (8, "6.1") => 8,
    (8, "7.1") => 680,
)

const TABLE_1_CATALOGUE = Catalogue(
    APNFunction(6, "1.1",
        monomial_expr(3)),
    APNFunction(6, "1.2",
        reference("1.1"),
        absolute_trace(trace_term(56, 3); scale = 1),
        relative_trace(3, trace_term(18, 9); scale = 1)),
    APNFunction(6, "2.3",
        monomial_expr(3), monomial_expr(1, 24; base = :u), monomial_expr(10)),
    APNFunction(6, "2.4",
        reference("2.1"),
        absolute_trace(trace_term(10, 3), trace_term(53, 5); scale = 3),
        relative_trace(3, trace_term(36, 9); scale = 3)),
    APNFunction(6, "2.7",
        reference("2.1"),
        absolute_trace(trace_term(34, 3), trace_term(48, 5)),
        relative_trace(3, trace_term(9, 9))),
    APNFunction(6, "2.10",
        reference("2.1"),
        absolute_trace(trace_term(24, 3), trace_term(28, 5); scale = 2),
        relative_trace(3, trace_term(0, 9); scale = 2)),
    APNFunction(6, "2.1",
        reference("2.3"),
        absolute_trace(trace_term(10, 3), trace_term(51, 5); scale = 42),
        relative_trace(3, trace_term(9, 9); scale = 42)),
    APNFunction(6, "2.5",
        reference("2.3"),
        absolute_trace(trace_term(31, 3), trace_term(49, 5); scale = 23),
        relative_trace(3, trace_term(9, 9); scale = 23)),
    APNFunction(6, "2.6",
        reference("2.3"),
        absolute_trace(trace_term(42, 3), trace_term(13, 5); scale = 12),
        relative_trace(3, trace_term(54, 9); scale = 12)),
    APNFunction(6, "2.8",
        reference("2.3"),
        absolute_trace(trace_term(51, 3), trace_term(60, 5); scale = 1),
        relative_trace(3, trace_term(18, 9); scale = 1)),
    APNFunction(6, "2.9",
        reference("2.3"),
        absolute_trace(trace_term(18, 3), trace_term(61, 5); scale = 14),
        relative_trace(3, trace_term(18, 9); scale = 14)),
    APNFunction(6, "2.11",
        reference("2.3"),
        absolute_trace(trace_term(50, 3), trace_term(56, 5); scale = 17)),
    APNFunction(6, "2.12",
        reference("2.3"),
        absolute_trace(trace_term(11, 3), trace_term(7, 5), trace_term(38, 7), trace_term(61, 11), trace_term(23, 13); scale = 19),
        relative_trace(3, trace_term(54, 9); scale = 19),
        relative_trace(2, trace_term(42, 21); scale = 19)),
    APNFunction(6, "2.2",
        reference("2.4"),
        absolute_trace(trace_term(54, 3), trace_term(47, 5); scale = 1),
        relative_trace(3, trace_term(9, 9); scale = 1)),

    APNFunction(8, "1.1",
        monomial_expr(3)),
    APNFunction(8, "1.2",
        reference("1.1"), absolute_trace(trace_term(48, 3), trace_term(0, 9))),
    APNFunction(8, "1.3",
        reference("1.1"), absolute_trace(trace_term(63, 3), trace_term(252, 9); scale = 1)),
    APNFunction(8, "1.4",
        reference("1.2"), absolute_trace(trace_term(84, 3), trace_term(213, 9); scale = 38)),
    APNFunction(8, "1.5",
        reference("1.2"), absolute_trace(trace_term(253, 3), trace_term(102, 9); scale = 51)),
    APNFunction(8, "1.6",
        reference("1.3"), absolute_trace(trace_term(68, 3), trace_term(235, 9); scale = 154)),
    APNFunction(8, "1.7",
        reference("1.4"), absolute_trace(trace_term(147, 3), trace_term(20, 9); scale = 69)),
    APNFunction(8, "1.8",
        reference("1.5"), absolute_trace(trace_term(153, 3), trace_term(51, 9); scale = 68)),
    APNFunction(8, "1.9",
        reference("1.6"), absolute_trace(trace_term(216, 3), trace_term(116, 9); scale = 35)),
    APNFunction(8, "1.10",
        reference("1.7"), absolute_trace(trace_term(232, 3), trace_term(195, 9); scale = 22)),
    APNFunction(8, "1.11",
        reference("1.8"), absolute_trace(trace_term(243, 3), trace_term(170, 9); scale = 85)),
    APNFunction(8, "1.12",
        reference("1.9"), absolute_trace(trace_term(172, 3), trace_term(31, 9); scale = 103)),
    APNFunction(8, "1.13",
        reference("1.10"), absolute_trace(trace_term(87, 3), trace_term(141, 5), trace_term(20, 9); scale = 90),
        relative_trace(4, trace_term(51, 17); scale = 90)),
    APNFunction(8, "1.14",
        reference("1.11"), absolute_trace(trace_term(160, 3), trace_term(250, 9); scale = 5)),
    APNFunction(8, "1.15",
        reference("1.11"), absolute_trace(trace_term(6, 3), trace_term(119, 9); scale = 102)),
    APNFunction(8, "1.16",
        reference("1.14"), absolute_trace(trace_term(133, 3), trace_term(30, 9); scale = 64)),
    APNFunction(8, "1.17",
        reference("1.16"), absolute_trace(trace_term(235, 3), trace_term(146, 9); scale = 78)),
    APNFunction(8, "2.1",
        monomial_expr(3), monomial_expr(17), monomial_expr(16, 18; base = :u), monomial_expr(16, 33; base = :u), monomial_expr(15, 48; base = :u)),
    APNFunction(8, "3.1",
        monomial_expr(3), monomial_expr(24, 6; base = :u), monomial_expr(182, 132; base = :u), monomial_expr(67, 192; base = :u)),
    APNFunction(8, "4.1",
        monomial_expr(3), monomial_expr(6), monomial_expr(68), monomial_expr(80), monomial_expr(132), monomial_expr(160)),
    APNFunction(8, "5.1",
        monomial_expr(3), monomial_expr(5), monomial_expr(18), monomial_expr(40), monomial_expr(66)),
    APNFunction(8, "6.1",
        monomial_expr(3), monomial_expr(12), monomial_expr(40), monomial_expr(66), monomial_expr(130)),
    APNFunction(8, "7.1",
        monomial_expr(57)),
)

const TABLE_1_FUNCTIONS = TABLE_1_CATALOGUE.functions

const TABLE_1_FUNCTION_BY_KEY = Dict((function_.n, function_.id) => function_
                                             for function_ in TABLE_1_FUNCTIONS)

const TABLE_1_EQUATION_IDS = Dict{Tuple{Int, String}, String}(
    (6, "1.1") => "1.1",
    (6, "1.2") => "1.2",
    (6, "2.3") => "2.1",
    (6, "2.4") => "2.2",
    (6, "2.7") => "2.3",
    (6, "2.10") => "2.4",
    (6, "2.1") => "2.5",
    (6, "2.5") => "2.6",
    (6, "2.6") => "2.7",
    (6, "2.8") => "2.8",
    (6, "2.9") => "2.9",
    (6, "2.11") => "2.10",
    (6, "2.12") => "2.11",
    (6, "2.2") => "2.12",
)

const TABLE_1_FUNCTION_BY_EQUATION_KEY = Dict((function_.n, get(TABLE_1_EQUATION_IDS, (function_.n, function_.id), function_.id)) => function_
                                                      for function_ in TABLE_1_FUNCTIONS)

const TABLE_1_GENERATION_CASES = [(n = function_.n, id = function_.id)
                                          for function_ in TABLE_1_FUNCTIONS]

const TABLE_1_CASES = [
    (n = 6, id = "1.1"),
    (n = 6, id = "1.2"),
    (n = 6, id = "2.1"),
    (n = 6, id = "2.2"),
    (n = 6, id = "2.3"),
    (n = 6, id = "2.4"),
    (n = 6, id = "2.5"),
    (n = 6, id = "2.6"),
    (n = 6, id = "2.7"),
    (n = 6, id = "2.8"),
    (n = 6, id = "2.9"),
    (n = 6, id = "2.10"),
    (n = 6, id = "2.11"),
    (n = 6, id = "2.12"),
    (n = 8, id = "1.1"),
    (n = 8, id = "1.2"),
    (n = 8, id = "1.3"),
    (n = 8, id = "1.4"),
    (n = 8, id = "1.5"),
    (n = 8, id = "1.6"),
    (n = 8, id = "1.7"),
    (n = 8, id = "1.8"),
    (n = 8, id = "1.9"),
    (n = 8, id = "1.10"),
    (n = 8, id = "1.11"),
    (n = 8, id = "1.12"),
    (n = 8, id = "1.13"),
    (n = 8, id = "1.14"),
    (n = 8, id = "1.15"),
    (n = 8, id = "1.16"),
    (n = 8, id = "1.17"),
    (n = 8, id = "2.1"),
    (n = 8, id = "3.1"),
    (n = 8, id = "4.1"),
    (n = 8, id = "5.1"),
    (n = 8, id = "6.1"),
    (n = 8, id = "7.1"),
]

function table_1_expected_permutations(n::Int, id::AbstractString)
    return TABLE_1_EXPECTED[(n, String(id))]
end

function table_1_formula(n::Int, id::AbstractString)
    return string(TABLE_1_FUNCTION_BY_KEY[(n, String(id))])
end

function table_1_equation_id(n::Int, id::AbstractString)
    key = (n, String(id))
    return get(TABLE_1_EQUATION_IDS, key, key[2])
end

function kaleyski_coefficient(field, exponent::Int)
    iszero(exponent) && return one(field)
    return gen(field)^exponent
end

function table_1_trace_value(field, x, trace_degree::Int, terms::Vector{APNTraceMonomial})
    value = zero(field)

    for trace_term in terms
        coefficient = kaleyski_coefficient(field, trace_term.coefficient_power)
        for i in 0:(trace_degree - 1)
            value += coefficient^(2^i) * x^(trace_term.exponent * 2^i)
        end
    end

    return value
end

function table_1_term_value(term::APNTerm, field, x)
    if term.coefficient isa OneCoefficient
        return x^term.exponent
    end

    if term.coefficient isa PowerCoefficient
        return kaleyski_coefficient(field, term.coefficient.exponent) * x^term.exponent
    end

    error("unsupported APN term coefficient: $(term.coefficient)")
end

function table_1_component_value(component::APNMonomial, field, x, x_int::Int, n::Int, generated_luts)
    return kaleyski_coefficient(field, component.coefficient_power) * x^component.exponent
end

function table_1_component_value(component::APNAbsoluteTrace, field, x, x_int::Int, n::Int, generated_luts)
    return kaleyski_coefficient(field, component.scale_power) *
           table_1_trace_value(field, x, n, component.terms)
end

function table_1_component_value(component::APNRelativeTrace, field, x, x_int::Int, n::Int, generated_luts)
    return kaleyski_coefficient(field, component.scale_power) *
           table_1_trace_value(field, x, component.extension_degree, component.terms)
end

function table_1_component_value(component::APNReference, field, x, x_int::Int, n::Int, generated_luts)
    function_ = TABLE_1_FUNCTION_BY_EQUATION_KEY[(n, component.id)]
    key = (n, function_.id)
    haskey(generated_luts, key) || error("referenced LUT $key has not been generated yet")
    return int_to_field_element(generated_luts[key][x_int + 1], field, n)
end

function table_1_generated_lut(case, generated_luts)
    field = GF(2, case.n, "a")
    function_ = TABLE_1_FUNCTION_BY_KEY[(case.n, case.id)]
    lut = Vector{Int}(undef, 2^case.n)

    for x_int in 0:(2^case.n - 1)
        x = int_to_field_element(x_int, field, case.n)
        value = zero(field)

        for term in function_.terms
            value += table_1_term_value(term, field, x)
        end

        for component in function_.components
            value += table_1_component_value(component, field, x, x_int, case.n, generated_luts)
        end

        lut[x_int + 1] = field_element_to_int(value, case.n)
    end

    return lut
end

function table_1_generated_luts()
    generated_luts = Dict{Tuple{Int, String}, Vector{Int}}()

    for case in TABLE_1_GENERATION_CASES
        generated_luts[(case.n, case.id)] = table_1_generated_lut(case, generated_luts)
    end

    return generated_luts
end

const TABLE_1_GENERATED_LUTS = table_1_generated_luts()

function table_1_context(; catalogue::Catalogue = TABLE_1_CATALOGUE)
    return ExperimentContext(
        :table_1,
        catalogue;
        description = "Reproduction context for Table 1 of Kaleyski, Deciding EA-equivalence via invariants.",
        source = "Kaleyski, Deciding EA-equivalence via invariants; test functions from Edel and Pott's APN catalogue.",
        root_dir = @__DIR__,
        artifacts_dir = joinpath(@__DIR__, "artifacts"),
        fixtures = Dict(:luts => TABLE_1_GENERATED_LUTS),
        expected = Dict(:external_linear_permutation_count => TABLE_1_EXPECTED),
        metadata = Dict(
            :case_order => TABLE_1_CASES,
            :generated_cases => TABLE_1_GENERATION_CASES,
            :equation_ids => TABLE_1_EQUATION_IDS,
        ),
    )
end

const TABLE_1_CONTEXT = table_1_context()

function table_1_selected_context(cases)
    functions_by_key = Dict(catalogue_key(function_) => function_
                            for function_ in TABLE_1_FUNCTIONS)
    functions = [functions_by_key[(case.n, case.id)] for case in cases]
    return table_1_context(catalogue = Catalogue(functions))
end

function table_1_cases(; dimensions = nothing, ids = nothing)
    cases = TABLE_1_CASES

    if dimensions !== nothing
        dimensions = dimensions isa Integer ? [dimensions] : dimensions
        wanted_dimensions = Set(dimensions)
        cases = [case for case in cases if case.n in wanted_dimensions]
    end

    if ids !== nothing
        ids = ids isa AbstractString ? [ids] : ids
        wanted_ids = Set(String.(ids))
        cases = [case for case in cases if case.id in wanted_ids]
    end

    return cases
end

function table_1_row(context::ExperimentContext,
                     function_::APNFunction;
                     parallel::Bool = Threads.nthreads() > 1,
                     log_level::Union{Symbol, AbstractString} = :quiet)
    found = 0

    elapsed = @elapsed begin
        lut = fixture(context, :luts)[catalogue_key(function_)]
        found = length(algorithm_1(lut,
                                   lut,
                                   function_.n;
                                   parallel = parallel,
                                   log_level = log_level))
    end

    expected = expected_value(context, :external_linear_permutation_count, function_)

    return (
        n = function_.n,
        id = function_.id,
        equation = table_1_equation_id(function_.n, function_.id),
        time = elapsed,
        found = found,
        expected = expected,
        parallel = parallel,
        log_level = Symbol(log_level),
        ok = found == expected,
    )
end

function run_table_1(; cases = TABLE_1_CASES,
                     dimensions = nothing,
                     ids = nothing,
                     parallel::Bool = Threads.nthreads() > 1,
                     log_level::Union{Symbol, AbstractString} = :quiet)
    selected_cases = cases === TABLE_1_CASES ? table_1_cases(dimensions = dimensions, ids = ids) : cases
    context = table_1_selected_context(selected_cases)
    rows = Any[]

    for function_ in context.catalogue.functions
        push!(rows, table_1_row(context, function_; parallel = parallel, log_level = log_level))
    end

    return ExperimentResult(
        ExperimentSpec(
            :external_linear_permutation_count,
            (context, function_) -> table_1_row(context, function_; parallel = parallel, log_level = log_level);
            description = "Number of external linear permutations reconstructed for each Table 1 function.",
            columns = [:n, :id, :equation, :time, :found, :expected, :delta, :parallel, :log_level, :ok],
        ),
        context.id,
        rows,
        RunMetadata(),
    )
end

function main_table_1(; kwargs...)
    result = run_table_1(; kwargs...)
    pretty_print_experiment_result(result;
                                   headers = Dict(:n => "n",
                                                  :id => "ID",
                                                  :log_level => "Log",
                                                  :ok => "OK"),
                                   alignments = Dict(:n => :right,
                                                     :time => :right,
                                                     :found => :right,
                                                     :expected => :right,
                                                     :delta => :right,
                                                     :parallel => :center,
                                                     :log_level => :center,
                                                     :ok => :center),
                                   mismatch_columns = [:n, :id, :equation, :found, :expected, :delta])
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_table_1()
end
