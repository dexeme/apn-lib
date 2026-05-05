julia --project=.

using Pkg
using Nemo
using APNLib
using JuliaInterpreter

julia > T = gen_permutation_tuple(n)
julia > generate_tuples_file(T, "tuples/AllTuples$n.h")
