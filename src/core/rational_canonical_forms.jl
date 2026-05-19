using Nemo

function integer_partitions(n::Int)
    n >= 0 || error("n must be non-negative")

    result = Vector{Vector{Int}}()

    function backtrack(remaining::Int, max_part::Int, current::Vector{Int})
        if remaining == 0
            push!(result, copy(current))
            return
        end

        for part in min(remaining, max_part):-1:1
            push!(current, part)
            backtrack(remaining - part, part, current)
            pop!(current)
        end
    end

    backtrack(n, n, Int[])
    return result
end

function monic_invertible_polynomial_coefficients(F, d::Int)
    d >= 1 || error("Dimension d must be positive")

    coeffs_list = Vector{Vector{typeof(one(F))}}()

    for value in 0:(2^(d + 1) - 1)
        bits = digits(value, base = 2, pad = d + 1)

        constant_term_is_one = bits[1] == 1
        leading_term_is_one = bits[d + 1] == 1

        if constant_term_is_one && leading_term_is_one
            push!(coeffs_list, [F(bit) for bit in bits])
        end
    end

    return coeffs_list
end

function companion_matrix_gf2(F, coeffs)
    d = length(coeffs) - 1

    d >= 1 || error("Degree must be at least 1")
    coeffs[1] == one(F) || error("Constant coefficient must be 1")
    coeffs[d + 1] == one(F) || error("Leading coefficient must be 1")

    companion = zero_matrix(F, d, d)

    for i in 2:d
        companion[i, i - 1] = one(F)
    end

    # Match Sage's companion_matrix(Q(list(v))) over GF(2):
    # the last column stores [c0, c1, ..., c_{d-1}]^T for
    # f(X) = c0 + c1*X + ... + c_{d-1}*X^(d-1) + X^d.
    for i in 1:d
        companion[i, d] = coeffs[i]
    end

    return companion
end

function blocks_for_rcf(d::Int)
    F = gf2()
    blocks = FqMatrix[]

    for coeffs in monic_invertible_polynomial_coefficients(F, d)
        push!(blocks, companion_matrix_gf2(F, coeffs))
    end

    return blocks
end

function polynomial_to_parent(poly, R)
    coeffs = [coeff(poly, i) for i in 0:degree(poly)]
    return R(coeffs)
end

function polynomial_equal(a, b)
    R = parent(b)

    a_same_parent = polynomial_to_parent(a, R)
    b_same_parent = polynomial_to_parent(b, R)

    return a_same_parent == b_same_parent
end

function matrix_is_similar(A::FqMatrix, B::FqMatrix)
    check_compatible_pair(A, B)

    order_A = matrix_multiplicative_order(A)
    order_B = matrix_multiplicative_order(B)

    if order_A != order_B
        return false
    end

    F = base_ring(A)
    n = nrows(A)
    identity = identity_matrix(F, n)

    # Order 1: both are identity matrices.
    if order_A == 1
        return true
    end

    # In characteristic 2, matrices of order 2 are unipotent.
    # Similarity is determined by the rank of A - I.
    if order_A == 2
        return rank(A - identity) == rank(B - identity)
    end

    # For the odd prime-order matrices used in this algorithm,
    # the matrices are semisimple, so the characteristic polynomial
    # determines the similarity class.
    char_A = charpoly(A)
    char_B = charpoly(B)

    return polynomial_equal(char_A, char_B)
end

function polynomial_divides(a, b)
    R = parent(b)

    a_same_parent = polynomial_to_parent(a, R)
    b_same_parent = polynomial_to_parent(b, R)

    _, r = divrem(b_same_parent, a_same_parent)

    return iszero(r)
end

function valid_rcf_block_sequence(blocks)
    isempty(blocks) && return false

    previous_poly = minpoly(blocks[1])

    for block in blocks
        current_poly = minpoly(block)

        if !polynomial_divides(current_poly, previous_poly)
            return false
        end

        previous_poly = current_poly
    end

    return true
end

function get_rcfs(n::Int)::Vector{FqMatrix}
    n >= 1 || error("n must be positive")

    matrices = FqMatrix[]

    for partition in integer_partitions(n)
        block_options = [blocks_for_rcf(part) for part in partition]

        for blocks in Base.Iterators.product(block_options...)
            if valid_rcf_block_sequence(blocks)
                reversed_blocks = FqMatrix[block for block in reverse(blocks)]
                push!(matrices, block_diagonal_gf2(reversed_blocks))
            end
        end
    end

    return matrices
end
