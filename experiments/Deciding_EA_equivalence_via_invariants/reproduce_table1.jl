using Nemo
using APNLib

include(joinpath(@__DIR__, "KaleyskiExperiments.jl"))
using .KaleyskiExperiments

function run_kaleyski_table1_case(case)
    context = kaleyski_table1_selected_context([case])
    result = run_kaleyski_table1_experiment(context = context)
    return only(result.rows)
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
    cases = KaleyskiExperiments.KALEYSKI_TABLE1_CASES

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
    context = kaleyski_table1_selected_context(cases)
    @info "Running Kaleyski Table 1 cases" cases = length(cases) threads = Threads.nthreads()
    print_markdown_header()

    result = run_kaleyski_table1_experiment(context = context)
    for (index, row) in pairs(result.rows)
        @info "Finished Kaleyski Table 1 case" index total = length(result.rows) n = row.n id = row.id
        print_markdown_row(row)
    end

    print_mismatch_summary(result.rows)
end

main()
