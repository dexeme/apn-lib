using Nemo

function block_diagonal_gf2(blocks::Vector{FqMatrix})::FqMatrix
    isempty(blocks) && error("At least one block is required")

    typed_blocks = FqMatrix[block for block in blocks]
    return block_diagonal_matrix(typed_blocks)
end

@doc"""
    matrix_multiplicative_order(A::FqMatrix) -> Int

Return the multiplicative order of the square matrix `A`, i.e. the smallest
positive integer `k` such that `A^k` is the identity matrix.
"""
function matrix_multiplicative_order(A::FqMatrix)::Int
    check_square(A)

    F = base_ring(A)
    n = nrows(A)
    I = identity_matrix(F, n)
    power = I
    max_iterations = 10_000_000

    for k in 1:max_iterations
        power = power * A

        if power == I
            return k
        end
    end

    error("Could not determine multiplicative order within $max_iterations iterations")
end
