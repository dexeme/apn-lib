# Run Tests

julia --project=. -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")'

# Parallel External Reconstruction

`reconstruct_external_linear_maps` uses the parallel backtracking implementation automatically when Julia is started with more than one thread.

Start Julia with threads:

```bash
julia -t auto --project=.
```

or choose a fixed number:

```bash
julia -t 8 --project=.
```

Run the external reconstruction:

```julia
using APNLib
using Nemo

n = 6
K = GF(2, n, "a")
R, x = Nemo.polynomial_ring(K, "x")

F = univariate_to_lut(x^3, n)
G = F

L1_candidates = reconstruct_external_linear_maps(F, G, n)
```

Enable progress logs:

```julia
L1_candidates = reconstruct_external_linear_maps(F, G, n, log_level = :info)
```

Available log levels are:

```julia
:quiet  # default
:info   # phase timings and parallel progress
:debug  # one line per first-level backtracking branch
```

Force one mode explicitly:

```julia
L1_parallel = reconstruct_external_linear_maps(F, G, n, parallel = true)
L1_serial = reconstruct_external_linear_maps(F, G, n, parallel = false)
```

Run one Kaleyski Table 1 case with threads:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=8 APNLIB_KALEYSKI_TABLE1_IDS=1.2 julia -t auto --project=. experiments/kaleyski_table1/reproduce_table1.jl
```

Run all Kaleyski Table 1 cases for one dimension:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=6 julia --project=. experiments/kaleyski_table1/reproduce_table1.jl
```

Enable algorithm logs from the script:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=6 APNLIB_LOG_LEVEL=info julia --project=. experiments/kaleyski_table1/reproduce_table1.jl
```

Enable parallel execution from the script by starting Julia with threads:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=6 APNLIB_LOG_LEVEL=info julia -t auto --project=. experiments/kaleyski_table1/reproduce_table1.jl
```

# Generate Tuples

julia --project=.

using Pkg; using Nemo; using APNLib

julia > include("experiments/linearly_self_equivalent_apn_permutations/tuple_generation.jl")
julia > T = gen_permutation_tuples(n)
julia > generate_tuples_file(T, "experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuples$n.h")

# Extract matrices A, B from tuples 

A, B = precomputed_tuple_matrices(n,class)

# Apply proposition 4
proposition4_filter(A, B, n)

# Full APN Search Workflow

After generating `experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuples$n.h`, generate the Julia constants with the tuple LUT and the extracted matrices:

```julia
using APNLib

generate_tuple_matrix_constants_file(n)
```

This creates:

```text
experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuplesMatrices$n.jl
```

Each class block contains:

```julia
const ALL_TUPLES_{n}_{class}_TUPLE = Int[...]
const ALL_TUPLES_{n}_{class}_A = Int[...]
const ALL_TUPLES_{n}_{class}_B = Int[...]
```

Load one class and run the proposition filters:

```julia
using APNLib
include("experiments/linearly_self_equivalent_apn_permutations/tuple_generation.jl")

n = 7
class_index = 1

A, B = precomputed_tuple_matrices(n, class_index)

passes_proposition4 = proposition4_filter(A, B, n)
passes_proposition5 = proposition5_filter(A, B, n)
```

Run Algorithm 1 and save the first S-box found back into the same matrix constants file:

```julia
include("experiments/linearly_self_equivalent_apn_permutations/search.jl")

solutions = APNSearch(
    n,
    A,
    B,
    max_solutions = 1,
    on_solution = sbox -> println(sbox),
    save_results = true,
    class_index = class_index,
)
```

When `save_results = true`, the result is written to:

```text
experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuplesMatrices$n.jl
```

with this format:

```julia
const ALL_TUPLES_{n}_{class}_SEARCH = Int[...]
```

For example:

```julia
const ALL_TUPLES_7_1_SEARCH = Int[...]
```

Run Algorithm 1 for multiple classes:

```julia
using APNLib
include("experiments/linearly_self_equivalent_apn_permutations/search.jl")

n = 7

# One class, an explicit class list, or all precomputed classes.
c = "all"
# c = [1, 4, 7, 10]

results = APNSearchClasses(
    n,
    c,
    excluded_class_indices = [13, 18, 21],
    max_solutions = 1,
    save_results = true,
    on_solution = (class_index, sbox) -> println("class $class_index finished"),
)
```

`APNSearchClasses` returns a `Dict{Int, Vector{Vector{Int}}}` mapping each class to its solutions. With `save_results = true`, each solution is written as `ALL_TUPLES_{n}_{class}_SEARCH`.

# Print Formatted Polynomials

Load the saved S-box constant and print the interpolated polynomial:

```julia
include("experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuplesMatrices7.jl")

sbox = ALL_TUPLES_7_1_SEARCH
println(format_sbox_polynomial(sbox, 7))
```

You can also print the polynomial directly from a solution returned by the search:

```julia
sbox = solutions[1]
println(format_sbox_polynomial(sbox, n))
```
