using Test
using Nemo
using APNLib

include("test_basicUtilities.jl")
include("test_tuplesGen.jl")
include("test_APNSearch.jl")
include("test_proposition4.jl")
include("test_proposition5.jl")
include("test_algorithm1.jl")
include("test_rcfs.jl")
include("test_families.jl")
include("test_trace_multiplicities.jl")
include("test_representations.jl")
isfile(joinpath(@__DIR__, "test_polynomial_json.jl")) && include("test_polynomial_json.jl")
