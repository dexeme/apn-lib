# Algorithm 1 General framework for reconstructing the outer permutation.
# -----------------------------------------------------------------------
# Input : Two (n, m)-functions F and G
# Output: All linear permutations L1 of (F_2)^m respecting the partitions induced by F and G
#
# for s in (F_2)^m do
#     compute M^F_4(0, s)
#     compute M^G_4(0, s)
# end
# partition (F_2)^m by equal values of M^F_4(0, s)
# partition (F_2)^m by equal values of M^G_4(0, s)
# if the partitions cannot be aligned by multiplicity and block size, return empty
# choose a basis B and possible image values U
# return Algorithm 2(K, C, B, U)

using APNLib

struct KaleyskiMultiplicityPartition
    blocks::Vector{Vector{Int}}
    block_index::Vector{Int}
    multiplicities::Vector{Int}
end

struct KaleyskiExternalReconstructionData
    left::KaleyskiMultiplicityPartition
    right::KaleyskiMultiplicityPartition
    basis::Vector{Int}
    candidate_values::Vector{Int}
end

function kaleyski_lut_from_table(table::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    APNLib.check_lut_values(table, n)
    return Int.(table)
end

function kaleyski_lut_from_table(table::AbstractDict{<:Integer, <:Integer}, n::Int)::Vector{Int}
    field_size = APNLib.space_size(n)
    lut = fill(-1, field_size)

    for (x, y) in table
        APNLib.check_space_value(x, n, name = "input values")
        APNLib.check_space_value(y, n, name = "output values")
        lut[Int(x) + 1] == -1 || error("table contains a repeated input: $x")
        lut[Int(x) + 1] = Int(y)
    end

    all(value -> value != -1, lut) || error("table is missing at least one input")
    return lut
end

function kaleyski_partition_from_multiplicities(multiplicities_by_value, n::Int)::KaleyskiMultiplicityPartition
    field_size = APNLib.space_size(n)
    values_by_multiplicity = Dict{Int, Vector{Int}}()

    for value in 0:(field_size - 1)
        haskey(multiplicities_by_value, value) || error("multiplicities are missing value $value")
        multiplicity = Int(multiplicities_by_value[value])
        push!(get!(values_by_multiplicity, multiplicity, Int[]), value)
    end

    multiplicities = sort!(collect(keys(values_by_multiplicity)))
    blocks = Vector{Vector{Int}}(undef, length(multiplicities))
    block_index = zeros(Int, field_size)

    for (index, multiplicity) in pairs(multiplicities)
        block = sort!(values_by_multiplicity[multiplicity])
        blocks[index] = block

        for value in block
            block_index[value + 1] = index
        end
    end

    return KaleyskiMultiplicityPartition(blocks, block_index, multiplicities)
end

function kaleyski_partition_by_multiplicity(F, n::Int)::KaleyskiMultiplicityPartition
    lut = kaleyski_lut_from_table(F, n)
    multiplicities = APNLib.multiplicities_sigma(lut, n, 4)
    return kaleyski_partition_from_multiplicities(multiplicities, n)
end

function kaleyski_gf2_rank(values::AbstractVector{<:Integer}, n::Int)::Int
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

function kaleyski_gf2_add_pivot!(pivots::Vector{Int}, value::Int)::Bool
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

function kaleyski_gf2_is_independent(value::Int, pivots::Vector{Int})::Bool
    reduced = value

    while reduced != 0
        bit_index = 63 - leading_zeros(reduced)
        pivot = pivots[bit_index + 1]
        pivot == 0 && return true
        reduced = xor(reduced, pivot)
    end

    return false
end

function kaleyski_gf2_basis(values::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    pivots = zeros(Int, n)
    basis = Int[]

    for raw_value in values
        value = Int(raw_value)

        if kaleyski_gf2_add_pivot!(pivots, value)
            push!(basis, value)
            length(basis) == n && return basis
        end
    end

    return basis
end

function kaleyski_aligned_partitions(left::KaleyskiMultiplicityPartition,
                                     right::KaleyskiMultiplicityPartition)
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

    return KaleyskiMultiplicityPartition(blocks, block_index, copy(left.multiplicities))
end

function kaleyski_select_basis_and_candidates(partition::KaleyskiMultiplicityPartition,
                                              right::KaleyskiMultiplicityPartition,
                                              n::Int)
    ordered_indices = sortperm(partition.blocks, by = length)
    block_count = length(ordered_indices)
    suffix_values = [Int[] for _ in 1:(block_count + 1)]

    for position in block_count:-1:1
        suffix_values[position] = vcat(partition.blocks[ordered_indices[position]], suffix_values[position + 1])
    end

    best_size = typemax(Int)
    best_indices = Int[]

    function search(position::Int, chosen::Vector{Int}, union_values::Vector{Int})
        length(union_values) >= best_size && return

        if kaleyski_gf2_rank(union_values, n) == n
            best_size = length(union_values)
            empty!(best_indices)
            append!(best_indices, chosen)
            return
        end

        position > block_count && return
        kaleyski_gf2_rank(vcat(union_values, suffix_values[position]), n) == n || return

        block_index = ordered_indices[position]
        push!(chosen, block_index)
        search(position + 1, chosen, vcat(union_values, partition.blocks[block_index]))
        pop!(chosen)
        search(position + 1, chosen, union_values)
        return
    end

    search(1, Int[], Int[])
    isempty(best_indices) && error("no union of multiplicity blocks contains a GF(2)-basis")

    basis_values = sort!(reduce(vcat, (partition.blocks[index] for index in best_indices), init = Int[]))
    basis = kaleyski_gf2_basis(basis_values, n)
    length(basis) == n || error("selected union does not contain a full basis")

    candidate_values = sort!(reduce(vcat, (right.blocks[index] for index in best_indices), init = Int[]))
    return basis, candidate_values
end

function kaleyski_prepare_external_reconstruction(F, G, n::Int)
    left = kaleyski_partition_by_multiplicity(F, n)
    right = kaleyski_partition_by_multiplicity(G, n)
    aligned_right = kaleyski_aligned_partitions(left, right)
    aligned_right === nothing && return nothing

    basis, candidate_values = kaleyski_select_basis_and_candidates(left, aligned_right, n)
    return KaleyskiExternalReconstructionData(left, aligned_right, basis, candidate_values)
end

function algorithm_1(F, G, n::Int; parallel::Bool = Threads.nthreads() > 1, log_level = :quiet)
    data = kaleyski_prepare_external_reconstruction(F, G, n)
    data === nothing && return Vector{Vector{Int}}()

    return algorithm_2(data.left,
                       data.right,
                       data.basis,
                       data.candidate_values,
                       n;
                       parallel = parallel,
                       log_level = log_level)
end
