JULIA ?= julia
ROOT := ../..
JULIA_PROJECT := --project=$(ROOT)

DIMENSIONS ?=
IDS ?=
PARALLEL ?=
LOG_LEVEL ?= quiet

.PHONY: instantiate tests test table-1 reproduce

instantiate:
	$(JULIA) $(JULIA_PROJECT) -e 'using Pkg; Pkg.instantiate()'

tests: instantiate
	$(JULIA) $(JULIA_PROJECT) -e 'include("KaleyskiExperiments.jl"); using .KaleyskiExperiments; @assert length(TABLE_1_CASES) == length(TABLE_1_EXPECTED); context = table_1_selected_context(table_1_cases(dimensions = 6, ids = "1.1")); @assert length(context.catalogue.functions) == 1'

test: tests

table-1: instantiate
	APNLIB_TABLE1_DIMENSIONS="$(DIMENSIONS)" APNLIB_TABLE1_IDS="$(IDS)" APNLIB_TABLE1_PARALLEL="$(PARALLEL)" APNLIB_LOG_LEVEL="$(LOG_LEVEL)" $(JULIA) $(JULIA_PROJECT) -e 'include("KaleyskiExperiments.jl"); using .KaleyskiExperiments; dimensions_value = get(ENV, "APNLIB_TABLE1_DIMENSIONS", ""); ids_value = get(ENV, "APNLIB_TABLE1_IDS", ""); parallel_value = get(ENV, "APNLIB_TABLE1_PARALLEL", ""); dimensions = isempty(dimensions_value) ? nothing : parse.(Int, split(dimensions_value, ",")); ids = isempty(ids_value) ? nothing : split(ids_value, ","); parallel = isempty(parallel_value) ? Base.Threads.nthreads() > 1 : parse(Bool, parallel_value); main_table_1(dimensions = dimensions, ids = ids, parallel = parallel, log_level = Symbol(get(ENV, "APNLIB_LOG_LEVEL", "quiet")))'

reproduce: table-1
