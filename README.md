# Run Tests

julia --project=. -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")'

# Generate Tuples

julia --project=.

using Pkg
using Nemo
using APNLib
using JuliaInterpreter

julia > T = gen_permutation_tuple(n)
julia > generate_tuples_file(T, "tuples/AllTuples$n.h")

# Extract matrices A, B from tuples 

A, B = precomputed_tuple_matrices(n,class)

# Apply proposition 4
proposition4_filter(A, B, n)

# Full APN Search Workflow

After generating `tuples/AllTuples$n.h`, generate the Julia constants with the tuple LUT and the extracted matrices:

```julia
using APNLib

generate_tuple_matrix_constants_file(n)
```

This creates:

```text
tuples/AllTuplesMatrices$n.jl
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

n = 7
class_index = 1

A, B = precomputed_tuple_matrices(n, class_index)

passes_proposition4 = proposition4_filter(A, B, n)
passes_proposition5 = proposition5_filter(A, B, n)
```

Run Algorithm 1 and save the first S-box found back into the same matrix constants file:

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

When `save_results = true`, the result is written to:

```text
tuples/AllTuplesMatrices$n.jl
```

with this format:

```julia
const ALL_TUPLES_{n}_{class}_SEARCH = Int[...]
```

For example:

```julia
const ALL_TUPLES_7_1_SEARCH = Int[...]
```

# Print Formatted Polynomials

Load the saved S-box constant and print the interpolated polynomial:

```julia
include("tuples/AllTuplesMatrices7.jl")

sbox = ALL_TUPLES_7_1_SEARCH
println(format_sbox_polynomial(sbox, 7))
```

You can also print the polynomial directly from a solution returned by the search:

```julia
sbox = solutions[1]
println(format_sbox_polynomial(sbox, n))
```

For the monomial representative check script, use:

```julia
include("src/apagar/teste.jl")

find_monomial_representative_for_class(1, 7)
```

This prints a representative in the paper format:

```text
x |-> x^5
```
