using Nemo

function check_square(A::FqMatrix; name::String = "Matrix")::Bool
    nrows(A) == ncols(A) || error("$name must be square")
end

function check_same_size(A::FqMatrix, B::FqMatrix)::Bool
    nrows(A) == nrows(B) || error("Matrices must have the same number of rows")
    ncols(A) == ncols(B) || error("Matrices must have the same number of columns")
end

function check_same_field(A::FqMatrix, B::FqMatrix)
    base_ring(A) == base_ring(B) || error("Matrices must be over the same field")
end

function check_compatible_pair(A::FqMatrix, B::FqMatrix):Bool
    check_square(A, name = "A")
    check_square(B, name = "B")
    check_same_size(A, B)
    check_same_field(A, B)
end

function block_diagonal_gf2(blocks)
    isempty(blocks) && error("At least one block is required")

    typed_blocks = FqMatrix[block for block in blocks]

    return block_diagonal_matrix(typed_blocks)
end