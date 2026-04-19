# finite_field_vectors.jl
# Julia Script

function int_to_column_vector(F, value::Int, n::Int)
    bits = digits(value, base = 2, pad = n)
    return matrix(F, n, 1, bits)
end

function each_column_vector(F, n::Int)
    return (int_to_column_vector(F, value, n) for value in 0:(2^n - 1))
end

function gf2_element_to_int(x)
    return iszero(x) ? 0 : 1
end

function column_vector_to_int(v::FqMatrix)
    n = nrows(v)
    value = 0

    for i in 1:n
        value += gf2_element_to_int(v[i, 1]) * 2^(i - 1)
    end

    return value
end