JULIA ?= julia
JULIA_PROJECT := --project=.

.PHONY: instantiate tests test experiment-tests all-tests kaleyski-table1-tests linearly-self-equivalent-tests

instantiate:
	$(JULIA) $(JULIA_PROJECT) -e 'using Pkg; Pkg.instantiate()'

tests: instantiate
	$(JULIA) $(JULIA_PROJECT) -e 'include("test/runtests.jl")'

test: tests

experiment-tests: kaleyski-table1-tests linearly-self-equivalent-tests

all-tests: tests experiment-tests

kaleyski-table1-tests:
	$(MAKE) -C experiments/Deciding_EA_equivalence_via_invariants tests

linearly-self-equivalent-tests:
	$(MAKE) -C experiments/linearly_self_equivalent_apn_permutations tests
