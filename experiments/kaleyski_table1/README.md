# Kaleyski Table 1 Experiment

This experiment reproduces the implemented permutation counts from Table 1
of Kaleyski, "Deciding EA-equivalence via invariants", using LUT fixtures
generated from the formulas in `table1_data.jl`.

## Tests

From this directory:

```bash
make tests
```

From the repository root:

```bash
make -C experiments/kaleyski_table1 tests
```

The tests check that the generated LUT fixture exists for every configured
case and that `reconstruct_external_linear_maps` returns the expected count.

## Reproduce Table 1 Cases

Run all configured cases:

```bash
make reproduce
```

Run only selected dimensions:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=6 make reproduce
```

Run selected IDs:

```bash
APNLIB_KALEYSKI_TABLE1_DIMENSIONS=8 APNLIB_KALEYSKI_TABLE1_IDS=1.2 make reproduce
```

Enable logs:

```bash
APNLIB_LOG_LEVEL=info make reproduce
```

Use Julia threads:

```bash
JULIA='julia -t auto' APNLIB_KALEYSKI_TABLE1_DIMENSIONS=6 make reproduce
```

Parallel execution can be controlled with `APNLIB_PARALLEL=auto`, `true`, or
`false`.

## Regenerate LUT Fixtures

The fixture file `fixtures/table1_luts.jl` is generated from `table1_data.jl`.
Regenerate it with:

```bash
make generate-luts
```
