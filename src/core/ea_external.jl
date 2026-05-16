struct MultiplicityPartition
    blocks::Vector{Vector{Int}}
    block_index::Vector{Int}
    multiplicities::Vector{Int}
end

function lut_from_table(table::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    check_lut_values(table, n)
    return Int.(table)
end

function lut_from_table(table::AbstractDict{<:Integer, <:Integer}, n::Int)::Vector{Int}
    space_size = 2^n
    lut = fill(-1, space_size)

    for (x, y) in table
        0 <= x < space_size || error("input values must be between 0 and $(space_size - 1)")
        0 <= y < space_size || error("output values must be between 0 and $(space_size - 1)")
        lut[Int(x) + 1] == -1 || error("table contains a repeated input: $x")
        lut[Int(x) + 1] = Int(y)
    end

    all(value -> value != -1, lut) || error("table is missing at least one input")
    return lut
end

function partition_by_multiplicity(lut::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                   n::Int;
                                   k::Int = 4)::MultiplicityPartition
    normalized_lut = lut_from_table(lut, n)
    multiplicities_by_value = multiplicities_sigma(normalized_lut, n, k)
    values_by_multiplicity = Dict{Int, Vector{Int}}()

    for value in 0:(2^n - 1)
        multiplicity = multiplicities_by_value[value]
        push!(get!(values_by_multiplicity, multiplicity, Int[]), value)
    end

    multiplicities = sort!(collect(keys(values_by_multiplicity)))
    blocks = Vector{Vector{Int}}(undef, length(multiplicities))
    block_index = zeros(Int, 2^n)

    for (index, multiplicity) in pairs(multiplicities)
        block = sort!(values_by_multiplicity[multiplicity])
        blocks[index] = block

        for value in block
            block_index[value + 1] = index
        end
    end

    return MultiplicityPartition(blocks, block_index, multiplicities)
end

function gf2_rank(values::AbstractVector{<:Integer}, n::Int)::Int
    pivots = zeros(Int, n)
    rank_value = 0

    for raw_value in values
        value = Int(raw_value)

        while value != 0
            bit_index = 63 - leading_zeros(value)

            if pivots[bit_index + 1] == 0
                pivots[bit_index + 1] = value
                rank_value += 1
                break
            end

            value = xor(value, pivots[bit_index + 1])
        end
    end

    return rank_value
end

function gf2_is_independent(value::Int, pivots::Vector{Int})::Bool
    reduced = value

    while reduced != 0
        bit_index = 63 - leading_zeros(reduced)
        pivot = pivots[bit_index + 1]
        pivot == 0 && return true
        reduced = xor(reduced, pivot)
    end

    return false
end

function gf2_add_pivot!(pivots::Vector{Int}, value::Int)::Bool
    reduced = value

    while reduced != 0
        bit_index = 63 - leading_zeros(reduced)

        if pivots[bit_index + 1] == 0
            pivots[bit_index + 1] = reduced
            return true
        end

        reduced = xor(reduced, pivots[bit_index + 1])
    end

    return false
end

function gf2_basis(values::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    pivots = zeros(Int, n)
    basis = Int[]

    for raw_value in values
        value = Int(raw_value)

        if gf2_add_pivot!(pivots, value)
            push!(basis, value)
            length(basis) == n && return basis
        end
    end

    return basis
end

function select_minimal_basis_union(partition::MultiplicityPartition, n::Int)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}}
    ordered_indices = sortperm(partition.blocks, by = length)
    block_count = length(ordered_indices)
    suffix_values = [Int[] for _ in 1:(block_count + 1)]

    for position in block_count:-1:1
        suffix_values[position] = vcat(partition.blocks[ordered_indices[position]], suffix_values[position + 1])
    end

    best_size = typemax(Int)
    best_indices = Int[]

    function search(position::Int, chosen::Vector{Int}, union_values::Vector{Int})::Nothing
        current_size = length(union_values)
        current_size >= best_size && return

        if gf2_rank(union_values, n) == n
            best_size = current_size
            empty!(best_indices)
            append!(best_indices, chosen)
            return
        end

        position > block_count && return
        gf2_rank(vcat(union_values, suffix_values[position]), n) == n || return

        block_index = ordered_indices[position]
        block = partition.blocks[block_index]

        push!(chosen, block_index)
        search(position + 1, chosen, vcat(union_values, block))
        pop!(chosen)

        search(position + 1, chosen, union_values)
        return
    end

    search(1, Int[], Int[])
    isempty(best_indices) && error("no union of multiplicity blocks contains a GF(2)-basis")

    union_values = sort!(reduce(vcat, (partition.blocks[index] for index in best_indices), init = Int[]))
    basis = gf2_basis(union_values, n)
    length(basis) == n || error("selected union does not contain a full basis")

    return basis, union_values, sort!(copy(best_indices))
end

function aligned_partitions(left::MultiplicityPartition, right::MultiplicityPartition)::Union{Nothing, MultiplicityPartition}
    length(left.blocks) == length(right.blocks) || return nothing

    right_by_multiplicity = Dict{Int, Vector{Int}}()
    for (index, multiplicity) in pairs(right.multiplicities)
        right_by_multiplicity[multiplicity] = right.blocks[index]
    end

    blocks = Vector{Vector{Int}}(undef, length(left.blocks))
    block_index = zeros(Int, length(right.block_index))

    for (index, multiplicity) in pairs(left.multiplicities)
        haskey(right_by_multiplicity, multiplicity) || return nothing
        block = right_by_multiplicity[multiplicity]
        length(block) == length(left.blocks[index]) || return nothing
        blocks[index] = block

        for value in block
            block_index[value + 1] = index
        end
    end

    return MultiplicityPartition(blocks, block_index, copy(left.multiplicities))
end

function backtrack_external_linear_maps(left::MultiplicityPartition,
                                        right::MultiplicityPartition,
                                        basis::Vector{Int},
                                        candidate_values::Vector{Int},
                                        n::Int)::Vector{Vector{Int}}
    length(basis) == n || error("basis must contain n elements")
    space_size = 2^n
    image = fill(-1, space_size)
    image[1] = 0
    domain_span = [0]
    image_pivots = zeros(Int, n)
    results = Vector{Vector{Int}}()

    function assign_basis_image(level::Int)::Nothing
        if level > n
            push!(results, copy(image))
            return
        end

        basis_value = basis[level]

        for candidate in candidate_values
            gf2_is_independent(candidate, image_pivots) || continue

            valid = true
            new_domain_values = Vector{Int}(undef, length(domain_span))
            new_image_values = Vector{Int}(undef, length(domain_span))

            @inbounds for index in eachindex(domain_span)
                x = domain_span[index]
                y = xor(x, basis_value)
                y_image = xor(image[x + 1], candidate)

                if right.block_index[y_image + 1] != left.block_index[y + 1]
                    valid = false
                    break
                end

                new_domain_values[index] = y
                new_image_values[index] = y_image
            end

            valid || continue

            old_span_length = length(domain_span)
            old_pivots = copy(image_pivots)
            gf2_add_pivot!(image_pivots, candidate)

            @inbounds for index in eachindex(new_domain_values)
                image[new_domain_values[index] + 1] = new_image_values[index]
                push!(domain_span, new_domain_values[index])
            end

            assign_basis_image(level + 1)

            resize!(domain_span, old_span_length)
            image_pivots .= old_pivots

            @inbounds for value in new_domain_values
                image[value + 1] = -1
            end
        end

        return
    end

    assign_basis_image(1)
    return results
end

function reconstruct_external_linear_maps(F::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          G::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          n::Int;
                                          k::Int = 4)::Vector{Vector{Int}}
    left = partition_by_multiplicity(F, n, k = k)
    right = partition_by_multiplicity(G, n, k = k)
    aligned_right = aligned_partitions(left, right)
    aligned_right === nothing && return Vector{Vector{Int}}()

    basis, _, chosen_indices = select_minimal_basis_union(left, n)
    candidate_values = sort!(reduce(vcat, (aligned_right.blocks[index] for index in chosen_indices), init = Int[]))

    return backtrack_external_linear_maps(left, aligned_right, basis, candidate_values, n)
end
