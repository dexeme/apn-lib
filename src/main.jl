include("all_tuples.jl")

function generate_for_n(n::Int)
    println("Generating tuples for n = $n")

    elapsed = @elapsed begin
        tuples = gen_permutation_tuples(n)
        generate_tuples_file(tuples, "AllTuples$n.h")
    end

    println("Generated AllTuples$n.h with $(length(tuples)) tuples")
    println("Elapsed time: $(round(elapsed, digits = 2)) seconds")
    println()

    return tuples
end

T6 = generate_for_n(6)
T7 = generate_for_n(7)
T8 = generate_for_n(8)