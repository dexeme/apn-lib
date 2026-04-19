# extract_apn_matrix.jl
# Julia Script

function extract_apn_matrix(sbox::Vector{Int}, n::Int)
    basis = [2^i for i in 0:(n-1)]
    compact_matrix = Int[]

    f_0 = sbox[1]

    for i in 1:n
        for j in (i + 1):n
            x = basis[i]
            y = basis[j]

            a_ij = sbox[(x ⊻ y) + 1] ⊻ sbox[x + 1] ⊻ sbox[y + 1] ⊻ f_0

            push!(compact_matrix, a_ij)
        end
    end

    return compact_matrix
end