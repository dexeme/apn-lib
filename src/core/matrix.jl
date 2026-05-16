using Nemo

function space_size(n::Int)::Int
    n >= 0 || error("n must be non-negative")
    return 2^n
end

function check_length(values, expected::Int; name::String = "values", unit::String = "entries")::Bool
    length(values) == expected || error("$name must have $expected $unit")
    return true
end

function check_space_length(values, n::Int; name::String = "values", unit::String = "entries")::Bool
    return check_length(values, space_size(n), name = name, unit = unit)
end

function check_integer_range(value::Integer, minimum::Int, maximum::Int; name::String = "value")::Bool
    minimum <= value <= maximum || error("$name must be between $minimum and $maximum")
    return true
end

function check_space_value(value::Integer, n::Int; name::String = "value")::Bool
    return check_integer_range(value, 0, space_size(n) - 1, name = name)
end

function check_square(A::FqMatrix; name::String = "Matrix")::Bool
    nrows(A) == ncols(A) || error("$name must be square")
    return true
end

function check_square(A::AbstractMatrix; name::String = "Matrix")::Bool
    size(A, 1) == size(A, 2) || error("$name must be square")
    return true
end

function check_same_size(A::FqMatrix, B::FqMatrix)::Bool
    nrows(A) == nrows(B) || error("Matrices must have the same number of rows")
    ncols(A) == ncols(B) || error("Matrices must have the same number of columns")
    return true
end

function check_same_size(A::AbstractMatrix, B::AbstractMatrix)::Bool
    size(A) == size(B) || error("Matrices must have the same size")
    return true
end

function check_same_field(A::FqMatrix, B::FqMatrix)::Bool
    base_ring(A) == base_ring(B) || error("Matrices must be over the same field")
    return true
end

function check_n_by_n_matrix(A::FqMatrix, n::Int; name::String = "Matrix")::Bool
    check_square(A, name = name)
    nrows(A) == n || error("$name must have $n rows")
    ncols(A) == n || error("$name must have $n columns")
    return true
end

function check_n_by_n_matrix(A::AbstractMatrix, n::Int; name::String = "Matrix")::Bool
    check_square(A, name = name)
    size(A) == (n, n) || error("$name must be a $n x $n matrix")
    return true
end

function check_gf2_matrix(A::FqMatrix; name::String = "Matrix")::Bool
    base_ring(A) == gf2() || error("$name must be over GF(2)")
    return true
end

function check_gf2_matrix(A::FqMatrix, n::Int; name::String = "Matrix")::Bool
    check_n_by_n_matrix(A, n, name = name)
    check_gf2_matrix(A, name = name)
    return true
end

function check_compatible_pair(A::FqMatrix, B::FqMatrix)::Bool
    check_square(A, name = "A")
    check_square(B, name = "B")
    check_same_size(A, B)
    check_same_field(A, B)
    return true
end

function check_compatible_pair(A::AbstractMatrix, B::AbstractMatrix)::Bool
    check_square(A, name = "A")
    check_square(B, name = "B")
    check_same_size(A, B)
    return true
end

function check_sbox_space_size(sbox::AbstractVector{<:Integer}, n::Int; name::String = "sbox")::Bool
    return check_space_length(sbox, n, name = name)
end

function check_lut_values(lut::AbstractVector{<:Integer}, n::Int; name::String = "lut")::Bool
    check_space_length(lut, n, name = name)

    for value in lut
        check_space_value(value, n, name = "$name values")
    end

    return true
end

function check_sbox_ddt_sizes(sbox::AbstractVector{<:Integer}, ddt::Matrix{Int})::Bool
    expected_size = length(sbox)
    check_square(ddt, name = "ddt")
    size(ddt) == (expected_size, expected_size) || error("ddt must be $expected_size x $expected_size")
    return true
end

function block_diagonal_gf2(blocks::Vector{FqMatrix})::FqMatrix
    isempty(blocks) && error("At least one block is required")

    typed_blocks = FqMatrix[block for block in blocks]
    return block_diagonal_matrix(typed_blocks)
end
