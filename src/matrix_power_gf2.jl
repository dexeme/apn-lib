# matrix_power_gf2.jl
# Julia Script

using Nemo;
include("utils.jl")

function matrix_power_gf2(A::Matrix, i::Int)
    is_quadratic(A) || error("Matrix must be square")

    result = Matrix{Int}(I, nrows, ncols)
    base = A .% 2

    for _ in 1:i
        result = (result * base) .% 2
    end

    return result
end