# check_order_space.jl
# Julia Script

using Nemo;
include("utils.jl")

function check_order_space(A::FqMatrix, i::Int)
    is_square(A) || error("Matrix must be square")

    F = base_ring(A)
    n = nrows(A)
    A_power = A^i

    count = 0

    for v in each_column_vector(F, n)
        if A_power * v == v
            count += 1
        end
    end

    return log(2, count)
end