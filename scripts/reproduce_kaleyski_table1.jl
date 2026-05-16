using Nemo
using APNLib

include(joinpath(@__DIR__, "..", "test", "kaleyski_table1_data.jl"))
include(joinpath(@__DIR__, "..", "test", "fixtures", "kaleyski_table1_luts.jl"))

function run_kaleyski_table1_case(case)
    obtained = 0
    expected = kaleyski_table1_expected_permutations(case.n, case.id)

    elapsed = @elapsed begin
        lut = KALEYSKI_TABLE1_GENERATED_LUTS[(case.n, case.id)]
        results = reconstruct_external_linear_maps(lut, lut, case.n)
        obtained = length(results)
    end

    return (
        n = case.n,
        id = case.id,
        equation = kaleyski_table1_equation_id(case.n, case.id),
        time = elapsed,
        permutations = obtained,
        expected = expected,
        ok = obtained == expected,
    )
end

function print_markdown_table(rows)
    println("| n | ID | Equation | Time | Permutations | Expected | OK |")
    println("|---:|:---|:---|---:|---:|---:|:---:|")

    for row in rows
        print_markdown_row(row)
    end
end

function print_markdown_header()
    println("| n | ID | Equation | Time | Permutations | Expected | OK |")
    println("|---:|:---|:---|---:|---:|---:|:---:|")
    flush(stdout)
end

function print_markdown_row(row)
    time = string(round(row.time; digits = 6))
    ok = row.ok ? "yes" : "no"
    println("| $(row.n) | $(row.id) | $(row.equation) | $time | $(row.permutations) | $(row.expected) | $ok |")
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
    @info "Running Kaleyski Table 1 cases" cases = length(cases)
    print_markdown_header()

    rows = []
    for (index, case) in pairs(cases)
        @info "Starting Kaleyski Table 1 case" index total = length(cases) n = case.n id = case.id
        row = run_kaleyski_table1_case(case)
        push!(rows, row)
        print_markdown_row(row)
    end
end

main()
