JULIA ?= julia
ROOT := ../..
JULIA_PROJECT := --project=$(ROOT)

.PHONY: instantiate tests test

instantiate:
	$(JULIA) $(JULIA_PROJECT) -e 'using Pkg; Pkg.instantiate()'

tests: instantiate
	$(JULIA) $(JULIA_PROJECT) runtests.jl

test: tests
