using Nemo
using APNLib

include(joinpath(@__DIR__, "table1_data.jl"))
include(joinpath(@__DIR__, "fixtures", "table1_luts.jl"))

function env_log_level()
    return Symbol(get(ENV, "APNLIB_LOG_LEVEL", "quiet"))
end

function env_parallel()
    value = lowercase(get(ENV, "APNLIB_PARALLEL", "auto"))
    value == "auto" && return Threads.nthreads() > 1
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("APNLIB_PARALLEL must be auto, true, or false")
end

function run_kaleyski_table1_case(case)
    obtained = 0
    expected = kaleyski_table1_expected_permutations(case.n, case.id)
    log_level = env_log_level()
    parallel = env_parallel()

    elapsed = @elapsed begin
        lut = KALEYSKI_TABLE1_GENERATED_LUTS[(case.n, case.id)]
        results = reconstruct_external_linear_maps(lut,
                                                   lut,
                                                   case.n,
                                                   parallel = parallel,
                                                   log_level = log_level)
        obtained = length(results)
    end

    return (
        n = case.n,
        id = case.id,
        equation = kaleyski_table1_equation_id(case.n, case.id),
        time = elapsed,
        found = obtained,
        expected = expected,
        parallel = parallel,
        log_level = log_level,
        ok = obtained == expected,
    )
end

function print_markdown_table(rows)
    println("| n | ID | Equation | Time | Found | Expected | Delta | Parallel | Log | OK |")
    println("|---:|:---|:---|---:|---:|---:|---:|:---:|:---:|:---:|")

    for row in rows
        print_markdown_row(row)
    end
end

function print_markdown_header()
    println("| n | ID | Equation | Time | Found | Expected | Delta | Parallel | Log | OK |")
    println("|---:|:---|:---|---:|---:|---:|---:|:---:|:---:|:---:|")
    flush(stdout)
end

function print_markdown_row(row)
    time = string(round(row.time; digits = 6))
    ok = row.ok ? "yes" : "no"
    delta = row.found - row.expected
    println("| $(row.n) | $(row.id) | $(row.equation) | $time | $(row.found) | $(row.expected) | $delta | $(row.parallel) | $(row.log_level) | $ok |")
    flush(stdout)
end

function print_mismatch_summary(rows)
    mismatches = [row for row in rows if !row.ok]
    isempty(mismatches) && return

    println()
    println("Mismatches:")
    println("| n | ID | Equation | Found | Expected | Delta |")
    println("|---:|:---|:---|---:|---:|---:|")

    for row in mismatches
        delta = row.found - row.expected
        println("| $(row.n) | $(row.id) | $(row.equation) | $(row.found) | $(row.expected) | $delta |")
    end

    flush(stdout)
end

function selected_cases()
    cases = KALEYSKI_TABLE1_CASES

    dimensions = get(ENV, "APNLIB_KALEYSKI_TABLE1_DIMENSIONS", "")
    if !isempty(dimensions)
        wanted = Set(parse.(Int, split(dimensions, ",")))
        cases = [case for case in cases if case.n in wanted]
    end

    ids = get(ENV, "APNLIB_KALEYSKI_TABLE1_IDS", "")
    if !isempty(ids)
        wanted = Set(strip.(split(ids, ",")))
        cases = [case for case in cases if case.id in wanted]
    end

    return cases
end

function main()
    cases = selected_cases()
    @info "Running Kaleyski Table 1 cases" cases = length(cases) threads = Threads.nthreads() parallel = env_parallel() log_level = env_log_level()
    print_markdown_header()

    rows = []
    for (index, case) in pairs(cases)
        @info "Starting Kaleyski Table 1 case" index total = length(cases) n = case.n id = case.id
        row = run_kaleyski_table1_case(case)
        push!(rows, row)
        print_markdown_row(row)
    end

    print_mismatch_summary(rows)
end

main()
