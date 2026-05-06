using Nemo

function _ensure_gf2_matrix(A::FqMatrix, n::Int; name::String = "Matrix")
    check_square(A, name = name)
    nrows(A) == n || error("$name must have $n rows")
    ncols(A) == n || error("$name must have $n columns")
    base_ring(A) == gf2() || error("$name must be over GF(2)")

    return A
end

function _ensure_gf2_matrix(A::AbstractMatrix, n::Int; name::String = "Matrix")
    size(A) == (n, n) || error("$name must be a $n x $n matrix")

    F = gf2()
    M = zero_matrix(F, n, n)

    for col in 1:n
        for row in 1:n
            M[row, col] = F(mod(Int(A[row, col]), 2))
        end
    end

    return M
end

function rank_gf2(M::FqMatrix)::Int
    base_ring(M) == gf2() || error("Matrix must be over GF(2)")

    return rank(M)
end

function nullity_gf2(M::FqMatrix)::Int
    base_ring(M) == gf2() || error("Matrix must be over GF(2)")

    return ncols(M) - rank_gf2(M)
end

function fixed_space_dimension_gf2(M_power::FqMatrix)::Int
    F = base_ring(M_power)
    n = nrows(M_power)
    I = identity_matrix(F, n)

    return nullity_gf2(M_power + I)
end

function filtro_proposicao_4(A::FqMatrix, B::FqMatrix, n::Int)::Bool
    A = _ensure_gf2_matrix(A, n, name = "A")
    B = _ensure_gf2_matrix(B, n, name = "B")
    check_compatible_pair(A, B)

    order_A = matrix_multiplicative_order(A)
    A_power = identity_matrix(base_ring(A), n)
    B_power = identity_matrix(base_ring(B), n)
    forbidden_dimensions = (2, 4, n - 1)

    for _ in 1:order_A
        A_power = A_power * A
        B_power = B_power * B

        d_A = fixed_space_dimension_gf2(A_power)
        d_B = fixed_space_dimension_gf2(B_power)

        if d_A != d_B || d_A in forbidden_dimensions || d_B in forbidden_dimensions
            return false
        end
    end

    return true
end

function filtro_proposicao_4(A::AbstractMatrix, B::AbstractMatrix, n::Int)::Bool
    A_gf2 = _ensure_gf2_matrix(A, n, name = "A")
    B_gf2 = _ensure_gf2_matrix(B, n, name = "B")

    return filtro_proposicao_4(A_gf2, B_gf2, n)
end

function check_order_space(A::FqMatrix, i::Int)
    check_square(A)

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

function matrix_to_sbox(A::FqMatrix)
    check_square(A)

    F = base_ring(A)
    n = nrows(A)

    table = Int[]

    for v in each_column_vector(F, n)
        push!(table, column_vector_to_int(A * v))
    end

    return table
end

function is_permutation_tuple(A::FqMatrix, B::FqMatrix)
    check_compatible_pair(A, B)

    n = nrows(A)

    for i in 0:(2^n - 1)
        if check_order_space(A, i) != check_order_space(B, i)
            return false
        end
    end

    return true
end

function same_multiplicative_order(A::FqMatrix, B::FqMatrix)
    return matrix_multiplicative_order(A) == matrix_multiplicative_order(B)
end

function has_prime_multiplicative_order(A::FqMatrix)::Bool
    order = matrix_multiplicative_order(A)
    return is_prime(order)
end

function is_power_similar(A_pair, B_pair)
    A1, A2 = A_pair
    B1, B2 = B_pair

    order_A1 = matrix_multiplicative_order(A1)
    order_A2 = matrix_multiplicative_order(A2)
    order_B1 = matrix_multiplicative_order(B1)
    order_B2 = matrix_multiplicative_order(B2)

    if order_A1 != order_B1 || order_A2 != order_B2
        return false
    end

    limit = max(order_A1, order_A2)

    for i in 0:(limit - 1)
        if matrix_is_similar(A1^i, B1) && matrix_is_similar(A2^i, B2)
            return true
        end
    end

    return false
end

function is_extended_power_similar(A_pair::Tuple{T,T}, B_pair::Tuple{T, T}) where T <: FqMatrix
    if is_power_similar(A_pair, B_pair)
        return true
    end

    A1, A2 = A_pair
    B1, B2 = B_pair

    return is_power_similar([inv(A1), inv(A2)], [B2, B1])
end

function candidate_pairs_with_same_order(matrices::Vector{T}) where T <: FqMatrix
    candidates = Vector{Tuple{T, T}}()
    for A in matrices
        for B in matrices
            if same_multiplicative_order(A, B)
                push!(candidates, (A, B))
            end
        end
    end
    return candidates
end


function remove_extended_power_similar_pairs(pairs::Vector{Tuple{T, T}}) where T <: FqMatrix
    n = length(pairs)
    active = trues(n)

    for i in 1:n
        !active[i] && continue
        for j in 1:n
            if i == j || !active[j]
                continue
            end
            if is_extended_power_similar(pairs[j], pairs[i])
                active[i] = false
                break
            end
        end
    end
    return pairs[active]
end

function filter_permutation_tuples(pairs::Vector{Tuple{T, T}}) where T <: FqMatrix
    valid_pairs = Tuple{T, T}[]

    for pair in pairs
        A, B = pair

        if is_permutation_tuple(A, B)
            push!(valid_pairs, pair)
        end
    end
    @show valid_pairs
    typeof(valid_pairs)
    return valid_pairs
end

# This script generates all possible matrix tuples for
# studying self-equivalence of permutations.
function gen_permutation_tuples(n::Int)::Vector{Tuple{FqMatrix, FqMatrix}}
    rcfs = get_rcfs(n)

    prime_order_rcfs = filter(has_prime_multiplicative_order, rcfs)

    candidates = candidate_pairs_with_same_order(prime_order_rcfs)

    representatives = remove_extended_power_similar_pairs(candidates)

    valid_pairs  = filter_permutation_tuples(representatives)

    return valid_pairs
end
