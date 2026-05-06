using Nemo

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

function block_diagonal_gf2(blocks::Vector{FqMatrix})::Vector{FqMatrix}
    isempty(blocks) && error("At least one block is required")

    typed_blocks = FqMatrix[block for block in blocks]
    return block_diagonal_matrix(typed_blocks)
end
