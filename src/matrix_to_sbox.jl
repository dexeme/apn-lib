# matrix_to_sbox.jl
# Julia Script

using Nemo;
include("utils.jl")

function matrix_to_sbox(A::FqMatrix)
    is_square(A) || error("Matrix must be square")

    F = base_ring(A)
    n = nrows(A)

    table = Int[]

    for v in each_column_vector(F, n)
        push!(table, column_vector_to_int(A * v))
    end

    return table
end