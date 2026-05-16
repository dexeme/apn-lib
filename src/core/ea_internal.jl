struct InternalReconstructionData
    domains::Vector{BitVector}
    basis::Vector{Int}
end

function bitset_values(set::BitVector)::Vector{Int}
    return [index - 1 for index in eachindex(set) if set[index]]
end

function build_o3_sets(lut::Vector{Int}, n::Int)::Vector{BitVector}
    field_size = space_size(n)
    sets = [falses(field_size) for _ in 1:field_size]

    @inbounds for x1 in 0:(field_size - 1)
        fx1 = lut[x1 + 1]

        for x2 in 0:(field_size - 1)
            x3 = xor(x1, x2)
            t = xor(fx1, lut[x2 + 1], lut[x3 + 1])
            set = sets[t + 1]
            set[x1 + 1] = true
            set[x2 + 1] = true
            set[x3 + 1] = true
        end
    end

    return sets
end

function restrict_internal_domains(F::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                   G::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                   n::Int)::Vector{BitVector}
    f_lut = lut_from_table(F, n)
    g_lut = lut_from_table(G, n)
    f_lut[1] == 0 || error("F(0) must be 0")
    g_lut[1] == 0 || error("G(0) must be 0")

    field_size = space_size(n)
    domains = [trues(field_size) for _ in 1:field_size]
    o3_sets = build_o3_sets(f_lut, n)

    @inbounds for y1 in 0:(field_size - 1)
        gy1 = g_lut[y1 + 1]

        for y2 in 0:(field_size - 1)
            y3 = xor(y1, y2)
            t = xor(gy1, g_lut[y2 + 1], g_lut[y3 + 1])
            allowed = o3_sets[t + 1]
            domains[y1 + 1] .&= allowed
            domains[y2 + 1] .&= allowed
            domains[y3 + 1] .&= allowed
        end
    end

    return domains
end

function optimized_internal_basis(domains::Vector{BitVector}, n::Int)::Vector{Int}
    field_size = space_size(n)
    check_length(domains, field_size, name = "domains", unit = "entries")

    ordered_values = sort!(collect(0:(field_size - 1)), by = value -> count(domains[value + 1]))
    pivots = zeros(Int, n)
    basis = Int[]

    for value in ordered_values
        value == 0 && continue

        if gf2_add_pivot!(pivots, value)
            push!(basis, value)
            length(basis) == n && return basis
        end
    end

    error("could not select a GF(2)-basis from the domain order")
end

function prepare_internal_reconstruction(F::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                         G::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                         n::Int)::InternalReconstructionData
    domains = restrict_internal_domains(F, G, n)
    basis = optimized_internal_basis(domains, n)
    return InternalReconstructionData(domains, basis)
end

function affine_lut(lut::Vector{Int}, n::Int)::Bool
    check_lut_values(lut, n)
    field_size = space_size(n)
    constant = lut[1]

    @inbounds for x in 0:(field_size - 1)
        ax = lut[x + 1]

        for y in 0:(field_size - 1)
            xor(ax, lut[y + 1], lut[xor(x, y) + 1], constant) == 0 || return false
        end
    end

    return true
end

function affine_on_span(linear_image::Vector{Int},
                        domain_span::Vector{Int},
                        F::Vector{Int},
                        G::Vector{Int},
                        c2::Int)::Bool
    a0 = xor(F[c2 + 1], G[1])

    @inbounds for x in domain_span
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

function reconstruct_internal_affine_maps(F::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          G::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          n::Int;
                                          is_quadratic::Bool = true,
                                          max_solutions::Union{Nothing, Int} = nothing)::Vector{Tuple{Vector{Int}, Vector{Int}}}
    f_lut = lut_from_table(F, n)
    g_lut = lut_from_table(G, n)
    data = prepare_internal_reconstruction(f_lut, g_lut, n)
    field_size = space_size(n)
    translations = is_quadratic ? [0] : collect(0:(field_size - 1))
    candidate_domains = [bitset_values(domain) for domain in data.domains]
    results = Vector{Tuple{Vector{Int}, Vector{Int}}}()

    function enough_results()::Bool
        return max_solutions !== nothing && length(results) >= max_solutions
    end

    for c2 in translations
        data.domains[1][c2 + 1] || continue

        linear_image = fill(-1, field_size)
        linear_image[1] = 0
        domain_span = [0]
        image_pivots = zeros(Int, n)

        function assign_basis_image(level::Int)::Nothing
            enough_results() && return

            if level > n
                a2_lut = Vector{Int}(undef, field_size)
                a_lut = Vector{Int}(undef, field_size)

                @inbounds for x in 0:(field_size - 1)
                    a2_x = xor(linear_image[x + 1], c2)
                    a2_lut[x + 1] = a2_x
                    a_lut[x + 1] = xor(f_lut[a2_x + 1], g_lut[x + 1])
                end

                if affine_lut(a_lut, n)
                    push!(results, (a2_lut, a_lut))
                end

                return
            end

            basis_value = data.basis[level]

            for candidate_a2 in candidate_domains[basis_value + 1]
                candidate_linear = xor(candidate_a2, c2)
                gf2_is_independent(candidate_linear, image_pivots) || continue

                valid = true
                old_span_length = length(domain_span)
                new_domain_values = Vector{Int}(undef, old_span_length)
                new_image_values = Vector{Int}(undef, old_span_length)

                @inbounds for index in 1:old_span_length
                    x = domain_span[index]
                    y = xor(x, basis_value)
                    y_linear = xor(linear_image[x + 1], candidate_linear)
                    y_a2 = xor(y_linear, c2)

                    if !data.domains[y + 1][y_a2 + 1]
                        valid = false
                        break
                    end

                    new_domain_values[index] = y
                    new_image_values[index] = y_linear
                end

                valid || continue

                old_pivots = copy(image_pivots)
                gf2_add_pivot!(image_pivots, candidate_linear)

                @inbounds for index in eachindex(new_domain_values)
                    linear_image[new_domain_values[index] + 1] = new_image_values[index]
                    push!(domain_span, new_domain_values[index])
                end

                if affine_on_span(linear_image, domain_span, f_lut, g_lut, c2)
                    assign_basis_image(level + 1)
                end

                resize!(domain_span, old_span_length)
                image_pivots .= old_pivots

                @inbounds for value in new_domain_values
                    linear_image[value + 1] = -1
                end
            end

            return
        end

        assign_basis_image(1)
        enough_results() && break
    end

    return results
end

algorithm3_reconstruct_internal = reconstruct_internal_affine_maps
