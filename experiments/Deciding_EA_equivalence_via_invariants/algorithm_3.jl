# Algorithm 3 Reconstructing the inner permutation A2.
# ----------------------------------------------------
# Input : Two (n, m)-functions F and G with F(0) = G(0) = 0
# Output: All affine permutations A2 such that F o A2 + G is affine.

using APNLib

function kaleyski_bitset_values(set::BitVector)::Vector{Int}
    return [index - 1 for index in eachindex(set) if set[index]]
end

function kaleyski_build_o3_sets(lut::Vector{Int}, n::Int)::Vector{BitVector}
    field_size = APNLib.space_size(n)
    sets = [falses(field_size) for _ in 1:field_size]

    for x1 in 0:(field_size - 1)
        for x2 in 0:(field_size - 1)
            x3 = xor(x1, x2)
            t = xor(lut[x1 + 1], lut[x2 + 1], lut[x3 + 1])
            sets[t + 1][x1 + 1] = true
            sets[t + 1][x2 + 1] = true
            sets[t + 1][x3 + 1] = true
        end
    end

    return sets
end

function kaleyski_restrict_internal_domains(F, G, n::Int)::Vector{BitVector}
    f_lut = kaleyski_lut_from_table(F, n)
    g_lut = kaleyski_lut_from_table(G, n)
    f_lut[1] == 0 || error("F(0) must be 0")
    g_lut[1] == 0 || error("G(0) must be 0")

    field_size = APNLib.space_size(n)
    domains = [trues(field_size) for _ in 1:field_size]
    o3_sets = kaleyski_build_o3_sets(f_lut, n)

    for x1 in 0:(field_size - 1)
        for x2 in 0:(field_size - 1)
            x3 = xor(x1, x2)
            t = xor(g_lut[x1 + 1], g_lut[x2 + 1], g_lut[x3 + 1])
            allowed = o3_sets[t + 1]
            domains[x1 + 1] .&= allowed
            domains[x2 + 1] .&= allowed
            domains[x3 + 1] .&= allowed
        end
    end

    return domains
end

function kaleyski_internal_basis(domains::Vector{BitVector}, n::Int)::Vector{Int}
    field_size = APNLib.space_size(n)
    ordered_values = sort!(collect(0:(field_size - 1)), by = value -> count(domains[value + 1]))
    pivots = zeros(Int, n)
    basis = Int[]

    for value in ordered_values
        value == 0 && continue

        if kaleyski_gf2_add_pivot!(pivots, value)
            push!(basis, value)
            length(basis) == n && return basis
        end
    end

    error("could not select a GF(2)-basis from the domain order")
end

function kaleyski_affine_lut(lut::Vector{Int}, n::Int)::Bool
    APNLib.check_lut_values(lut, n)
    field_size = APNLib.space_size(n)
    constant = lut[1]

    for x in 0:(field_size - 1)
        for y in 0:(field_size - 1)
            xor(lut[x + 1], lut[y + 1], lut[xor(x, y) + 1], constant) == 0 || return false
        end
    end

    return true
end

function kaleyski_affine_on_span(linear_image::Vector{Int},
                                 domain_span::Vector{Int},
                                 F::Vector{Int},
                                 G::Vector{Int},
                                 c2::Int)::Bool
    a0 = xor(F[c2 + 1], G[1])

    for x in domain_span
        ax = xor(F[xor(linear_image[x + 1], c2) + 1], G[x + 1])

        for y in domain_span
            xy = xor(x, y)
            ay = xor(F[xor(linear_image[y + 1], c2) + 1], G[y + 1])
            axy = xor(F[xor(linear_image[xy + 1], c2) + 1], G[xy + 1])
            xor(ax, ay, axy, a0) == 0 || return false
        end
    end

    return true
end

function algorithm_3(F, G, n::Int; is_quadratic::Bool = true, max_solutions = nothing)
    f_lut = kaleyski_lut_from_table(F, n)
    g_lut = kaleyski_lut_from_table(G, n)
    domains = kaleyski_restrict_internal_domains(f_lut, g_lut, n)
    basis = kaleyski_internal_basis(domains, n)
    field_size = APNLib.space_size(n)
    translations = is_quadratic ? [0] : collect(0:(field_size - 1))
    candidate_domains = [kaleyski_bitset_values(domain) for domain in domains]
    results = Vector{Tuple{Vector{Int}, Vector{Int}}}()

    for c2 in translations
        domains[1][c2 + 1] || continue

        L2 = fill(-1, field_size)
        L2[1] = 0
        domain_span = [0]
        image_pivots = zeros(Int, n)

        function assign(i::Int)::Bool
            if i > n
                A2 = Vector{Int}(undef, field_size)
                A = Vector{Int}(undef, field_size)

                for x in 0:(field_size - 1)
                    a2_x = xor(L2[x + 1], c2)
                    A2[x + 1] = a2_x
                    A[x + 1] = xor(f_lut[a2_x + 1], g_lut[x + 1])
                end

                if kaleyski_affine_lut(A, n)
                    push!(results, (A2, A))
                    return max_solutions === nothing || length(results) < max_solutions
                end

                return true
            end

            basis_value = basis[i]

            for candidate_a2 in candidate_domains[basis_value + 1]
                candidate_l2 = xor(candidate_a2, c2)
                kaleyski_gf2_is_independent(candidate_l2, image_pivots) || continue

                partition_preserved = true
                old_span_length = length(domain_span)
                new_domain_values = Vector{Int}(undef, old_span_length)
                new_image_values = Vector{Int}(undef, old_span_length)

                for index in 1:old_span_length
                    x = domain_span[index]
                    y = xor(x, basis_value)
                    y_l2 = xor(L2[x + 1], candidate_l2)
                    y_a2 = xor(y_l2, c2)

                    if !domains[y + 1][y_a2 + 1]
                        partition_preserved = false
                        break
                    end

                    new_domain_values[index] = y
                    new_image_values[index] = y_l2
                end

                partition_preserved || continue

                old_pivots = copy(image_pivots)
                kaleyski_gf2_add_pivot!(image_pivots, candidate_l2)

                for index in eachindex(new_domain_values)
                    L2[new_domain_values[index] + 1] = new_image_values[index]
                    push!(domain_span, new_domain_values[index])
                end

                should_continue = kaleyski_affine_on_span(L2, domain_span, f_lut, g_lut, c2) ? assign(i + 1) : true

                resize!(domain_span, old_span_length)
                image_pivots .= old_pivots

                for value in new_domain_values
                    L2[value + 1] = -1
                end

                should_continue || return false
            end

            return true
        end

        assign(1) || return results
    end

    return results
end
