JULIA ?= julia
ROOT := ../..
JULIA_PROJECT := --project=$(ROOT)

.PHONY: instantiate tests test reproduce generate-luts

instantiate:
	$(JULIA) $(JULIA_PROJECT) -e 'using Pkg; Pkg.instantiate()'

tests: instantiate
	$(JULIA) $(JULIA_PROJECT) runtests.jl

test: tests

reproduce: instantiate
	$(JULIA) $(JULIA_PROJECT) reproduce_table1.jl

generate-luts: instantiate
	$(JULIA) $(JULIA_PROJECT) generate_luts.jl
