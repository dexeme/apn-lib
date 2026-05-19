# Linearly Self-Equivalent APN Permutations Experiment

This experiment contains the tuple generation, proposition filters, and
Algorithm 1 search workflow for linearly self-equivalent APN permutations.

## Tests

From this directory:

```bash
make tests
```

From the repository root:

```bash
make -C experiments/linearly_self_equivalent_apn_permutations tests
```

The test suite loads the tuple-generation code, proposition filters, APN
search helpers, and Algorithm 1 regression tests.

## Generate Tuples

Start Julia from the repository root:

```bash
julia --project=.
```

Then load the experiment helpers:

```julia
using APNLib
using Nemo

include("experiments/linearly_self_equivalent_apn_permutations/tuple_generation.jl")
```

Generate permutation tuples for a dimension:

```julia
n = 6
T = gen_permutation_tuples(n)
```

Generate the Julia constants with tuple LUTs and extracted matrices:

```julia
generate_tuple_matrix_constants_file(n)
```

This writes:

```text
experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuplesMatrices$n.jl
```

## Proposition Filters

Load a precomputed class and apply the filters:

```julia
using APNLib
include("experiments/linearly_self_equivalent_apn_permutations/tuple_generation.jl")

n = 7
class_index = 1

A, B = precomputed_tuple_matrices(n, class_index)

passes_proposition4 = proposition4_filter(A, B, n)
passes_proposition5 = proposition5_filter(A, B, n)
```

## Algorithm 1 Search

Load the search implementation:

```julia
include("experiments/linearly_self_equivalent_apn_permutations/search.jl")
```

Run one class and save the first S-box found:

```julia
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

Run multiple classes:

```julia
results = APNSearchClasses(
    n,
    "all",
    excluded_class_indices = [13, 18, 21],
    max_solutions = 1,
    save_results = true,
    on_solution = (class_index, sbox) -> println("class $class_index finished"),
)
```

`APNSearchClasses` returns a `Dict{Int, Vector{Vector{Int}}}` mapping each
class to the solutions found.

## Print Saved Polynomials

Load a saved S-box constant and print its interpolated polynomial:

```julia
using APNLib

include("experiments/linearly_self_equivalent_apn_permutations/tuples/AllTuplesMatrices7.jl")

sbox = ALL_TUPLES_7_1_SEARCH
println(format_sbox_polynomial(sbox, 7))
```
