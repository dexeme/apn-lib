using Nemo

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
