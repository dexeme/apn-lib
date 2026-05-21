# Deciding EA-Equivalence via Invariants

This experiment implements the algorithms and the reproduction of Table 1 from
Nikolay Kaleyski's paper, "Deciding EA-equivalence via invariants".

The main files are:

- `KaleyskiExperiments.jl`: experiment module.
- `algorithm_1.jl`: reconstruction of the external linear permutation.
- `algorithm_2.jl`: search for linear permutations that respect partitions.
- `algorithm_3.jl`: reconstruction of the internal affine permutation.
- `table_1.jl`: data and routine for reproducing Table 1.

## Requirements

Use the repository's Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This directory's Makefile runs this step automatically before its targets.

## Makefile

From this directory:

```bash
make tests
make table-1
```

From the repository root:

```bash
make -C experiments/Deciding_EA_equivalence_via_invariants tests
make -C experiments/Deciding_EA_equivalence_via_invariants table-1
```

To run only some Table 1 entries:

```bash
make table-1 DIMENSIONS=6 IDS=1.1
make table-1 DIMENSIONS=8 IDS=1.1,1.2 LOG_LEVEL=info
```

It is also possible to control parallelization:

```bash
make table-1 DIMENSIONS=6 IDS=1.1 PARALLEL=false
julia -t auto --project=../.. -e 'include("KaleyskiExperiments.jl"); using .KaleyskiExperiments; main_table_1(dimensions = 6, ids = "1.1")'
```

## Julia Examples

Start Julia from the repository root:

```bash
julia --project=.
```

Load the experiment:

```julia
include("experiments/Deciding_EA_equivalence_via_invariants/KaleyskiExperiments.jl")
using .KaleyskiExperiments
```

List the available cases:

```julia
table_1_cases()
table_1_cases(dimensions = 6)
table_1_cases(dimensions = 8, ids = ["1.1", "1.2"])
```

Run a specific Table 1 entry:

```julia
result = main_table_1(dimensions = 6, ids = "1.1")
```

Run a subset with logs:

```julia
result = main_table_1(
    dimensions = 8,
    ids = ["1.1", "1.2"],
    parallel = true,
    log_level = :info,
)
```

Use Algorithm 1 directly on a LUT generated for a Table 1 case:

```julia
context = table_1_selected_context(table_1_cases(dimensions = 6, ids = "1.1"))
function_ = context.catalogue.functions[1]
lut = context.fixtures[:luts][(function_.n, function_.id)]

linear_maps = algorithm_1(lut, lut, function_.n)
length(linear_maps)
```

## Citation

```bibtex
@article{Kaleyski:2022,
  doi = {10.1007/s12095-021-00513-y},
  author = {Nikolay Kaleyski},
  title = {{Deciding {EA}-equivalence via invariants}},
  journal = {{Cryptography and Communications}},
  volume = {14},
  number = {2},
  pages = {271--290},
  year = {2022},
  month = mar,
}
```
