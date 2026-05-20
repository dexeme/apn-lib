# APNLib

APNLib is a Julia library for experiments with APN functions over binary
finite fields. It includes utilities for finite-field representations,
linear maps, EA-equivalence routines, APN lookup tables, polynomial
conversion, trace components, and regression experiments used in the TCC.

## Requirements

- Julia
- The Julia dependency declared in `Project.toml` (`Nemo`)

Install dependencies with:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The Makefiles run this step automatically before tests.

## Tests

Run the main library test suite from the repository root:

```bash
make tests
```

The same target is also available as:

```bash
make test
```

To run the experiment test suites from the root:

```bash
make experiment-tests
```

To run everything:

```bash
make all-tests
```

You can also run each experiment directly:

```bash
make -C experiments/kaleyski_table1 tests
make -C experiments/linearly_self_equivalent_apn_permutations tests
```

## Library Usage

Start Julia in the project:

```bash
julia --project=.
```

Load the package:

```julia
using APNLib
using Nemo
```

Example: reconstruct external linear maps for the Gold function over `GF(2^6)`:

```julia
n = 6

gold_polynomial = APNFunction(monomial_expr(3))
F = univariate_to_lut(gold_polynomial, n)
G = F

L1_candidates = reconstruct_external_linear_maps(F, G, n)
```

The external reconstruction routine can use Julia threads. Start Julia with:

```bash
julia -t auto --project=.
```

or with a fixed number of threads:

```bash
julia -t 8 --project=.
```

Progress logging is controlled by the `log_level` keyword:

```julia
L1_candidates = reconstruct_external_linear_maps(F, G, n, log_level = :info)
```

Available log levels are `:quiet`, `:info`, and `:debug`.

Parallel mode can also be selected explicitly:

```julia
L1_parallel = reconstruct_external_linear_maps(F, G, n, parallel = true)
L1_serial = reconstruct_external_linear_maps(F, G, n, parallel = false)
```

## Experiments

Experiment documentation is kept next to each experiment:

- [`experiments/kaleyski_table1/README.md`](experiments/kaleyski_table1/README.md)
- [`experiments/linearly_self_equivalent_apn_permutations/README.md`](experiments/linearly_self_equivalent_apn_permutations/README.md)
