using Nemo

@doc"""
    matrix_multiplicative_order(A::FqMatrix) -> Int
Returns the multiplicative order of the matrix `A`, i.e., the smallest positive integer `k` such that `A^k` is the identity matrix.

### Input
- `A::FqMatrix`: A square matrix over a finite field.

### Output
- `Int`: The multiplicative order of the matrix `A`.

@todo: optimize this function by using the characteristic polynomial and its factorization
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
