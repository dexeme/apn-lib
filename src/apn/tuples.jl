using Nemo

function _ensure_gf2_matrix(A::FqMatrix, n::Int; name::String = "Matrix")
    check_gf2_matrix(A, n, name = name)

    return A
end

function _ensure_gf2_matrix(A::AbstractMatrix, n::Int; name::String = "Matrix")
    check_n_by_n_matrix(A, n, name = name)

    F = gf2()
    matrix_gf2 = zero_matrix(F, n, n)

    for col in 1:n
        for row in 1:n
            matrix_gf2[row, col] = F(mod(Int(A[row, col]), 2))
        end
    end

    return matrix_gf2
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

function permutation_cycle_structure(permutation::Vector{Int})::Vector{Int}
    permutation_size = length(permutation)
    seen = falses(permutation_size)
    cycle_lengths = Int[]

    @inbounds for start in 0:(permutation_size - 1)
        seen[start + 1] && continue

        current = start
        cycle_length = 0

        while !seen[current + 1]
            seen[current + 1] = true
            cycle_length += 1
            current = permutation[current + 1]
            0 <= current < permutation_size || error("Permutation value out of range: $current")
        end

        push!(cycle_lengths, cycle_length)
    end

    return sort!(cycle_lengths)
end

function matrix_cycle_structure(A::FqMatrix)::Vector{Int}
    return permutation_cycle_structure(matrix_to_sbox(A))
end

function same_cycle_structure(A::FqMatrix, B::FqMatrix)::Bool
    return matrix_cycle_structure(A) == matrix_cycle_structure(B)
end
