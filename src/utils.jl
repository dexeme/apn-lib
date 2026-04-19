# int_to_column_vector.jl
# Julia Script

function is_quadratic(A::FqMatrix)
    return nrows(A) == ncols(A)
end

function fq2_elem_to_int(x)
    return iszero(x) ? 0 : 1
end

function int_to_column_vector(F, value::Int, n::Int)
    bits = digits(value, base = 2, pad = n)
    return matrix(F, n, 1, bits)
end

function column_vector_to_int(v::FqMatrix)
    n = nrows(v)

    value = 0

    for j in 1:n
        bit = fq2_elem_to_int(v[j, 1])
        value += bit * 2^(j - 1)
    end

    return value
end

function all_column_vectors(F, n::Int)
    return [int_to_column_vector(F, value, n) for value in 0:(2^n - 1)]
end