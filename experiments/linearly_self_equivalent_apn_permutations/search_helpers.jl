using Nemo
using APNLib

@doc"""
    even_hamming_weight_differences(space_size::Int) -> Vector{Int}

Return the nonzero differences whose binary representation has even Hamming
weight. The linearly self-equivalent APN search uses this restricted set when
it updates the partial difference distribution table.
"""
function even_hamming_weight_differences(space_size::Int)::Vector{Int}
    return [alpha for alpha in 1:(space_size - 1) if iseven(count_ones(alpha))]
end

function updateDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int}, delta::Int)::Bool
    delta in (-2, 2) || error("Unsupported DDT delta: $delta")

    c_value = sbox[c + 1]
    c_value != -1 || error("sbox[$c] must be assigned before updating the DDT")
    field_size = length(sbox)

    @inbounds for alpha in even_hamming_weight_differences(field_size)
        paired_x = xor(c, alpha)
        paired_value = sbox[paired_x + 1]
        paired_value == -1 && continue

        out_diff = xor(c_value, paired_value)
        ddt_value = ddt[alpha + 1, out_diff + 1] + delta
        ddt[alpha + 1, out_diff + 1] = ddt_value
        delta == 2 && ddt_value > 2 && return false
        delta == -2 && ddt_value == 2 && break
    end

    return true
end

@doc"""
    addDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int}) -> Bool

Add the contribution of the already assigned input `c` to a partial DDT.
Unassigned S-box entries must be represented by `-1`. The function returns
`false` as soon as an APN violation is detected, i.e. a DDT entry exceeds `2`.
"""
function addDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return addDDTInformationUnchecked(c, sbox, ddt)
end

function addDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    return updateDDTInformationUnchecked(c, sbox, ddt, 2)
end

@doc"""
    removeDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int}) -> Bool

Remove the contribution of the assigned input `c` from a partial DDT. This is
the inverse operation of [`addDDTInformation`](@ref) and is intended for
backtracking searches over partially assigned S-boxes.
"""
function removeDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return removeDDTInformationUnchecked(c, sbox, ddt)
end

function removeDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    return updateDDTInformationUnchecked(c, sbox, ddt, -2)
end

@doc"""
    isComplete(sbox::Vector{Int}) -> Bool

Return whether every entry of a partial S-box has been assigned. Unassigned
entries are represented by `-1`.
"""
function isComplete(sbox::Vector{Int})::Bool
    @inbounds for value in sbox
        value == -1 && return false
    end

    return true
end

@doc"""
    nextFreePosition(sbox::Vector{Int}) -> Int
    nextFreePosition(sbox::Vector{Int}, visit_order::Vector{Int}) -> Int

Return the next unassigned input of a partial S-box, or `-1` when the S-box is
complete. The two-argument form respects a caller-provided zero-based visit
order.
"""
function nextFreePosition(sbox::Vector{Int})::Int
    @inbounds for index in eachindex(sbox)
        sbox[index] == -1 && return index - 1
    end

    return -1
end

function nextFreePosition(sbox::Vector{Int}, visit_order::Vector{Int})::Int
    @inbounds for value in visit_order
        sbox[value + 1] == -1 && return value
    end

    return -1
end

@doc"""
    standard_visit_order(n::Int) -> Vector{Int}

Return the natural zero-based visit order `0:(2^n - 1)` for S-box searches.
"""
function standard_visit_order(n::Int)::Vector{Int}
    return collect(0:(space_size(n) - 1))
end

@doc"""
    offset_visit_order(n::Int, offset::Int) -> Vector{Int}

Return a cyclic shift of the natural zero-based visit order in `GF(2)^n`.
"""
function offset_visit_order(n::Int, offset::Int)::Vector{Int}
    field_size = space_size(n)
    return [mod(index + offset, field_size) for index in 0:(field_size - 1)]
end

@doc"""
    c_reference_visit_order(n::Int, class_index::Union{Nothing, Int}=nothing) -> Vector{Int}

Return the visit order used to reproduce the C reference search. Most classes
use the natural order; the known exceptional `n = 7` classes use fixed offsets.
"""
function c_reference_visit_order(n::Int, class_index::Union{Nothing, Int} = nothing)::Vector{Int}
    if n == 7 && class_index == 24
        return offset_visit_order(n, 16)
    elseif n == 7 && class_index == 26
        return offset_visit_order(n, 8)
    end

    return standard_visit_order(n)
end

@doc"""
    int_matrix_to_lut(M::AbstractMatrix, n::Int) -> Vector{Int}

Convert an `n x n` integer matrix, interpreted modulo two, into the lookup
table of the corresponding linear map on `GF(2)^n`.
"""
function int_matrix_to_lut(M::AbstractMatrix, n::Int)::Vector{Int}
    check_n_by_n_matrix(M, n)
    field_size = space_size(n)
    table = Vector{Int}(undef, field_size)

    @inbounds for x in 0:(field_size - 1)
        output = 0

        for row in 1:n
            output_bit = false

            for col in 1:n
                if isodd((x >> (col - 1)) & 1) && isodd(Int(M[row, col]))
                    output_bit = !output_bit
                end
            end

            if output_bit
                output |= 1 << (row - 1)
            end
        end

        table[x + 1] = output
    end

    return table
end

@doc"""
    linear_map_lut(M, n::Int) -> Vector{Int}

Convert a binary `n x n` matrix into the lookup table of its linear map on
zero-based vector encodings. Nemo matrices must be over `GF(2)`; integer
matrices are reduced modulo two.
"""
function linear_map_lut(M::FqMatrix, n::Int)::Vector{Int}
    check_gf2_matrix(M, n)
    return matrix_to_sbox(M)
end

function linear_map_lut(M::AbstractMatrix, n::Int)::Vector{Int}
    return int_matrix_to_lut(M, n)
end

@doc"""
    orbit_orders(apply_map::Vector{Int}) -> Vector{Int}

Return a vector whose entry `x + 1` is the length of the orbit of `x` under the
zero-based map table `apply_map`.
"""
function orbit_orders(apply_map::Vector{Int})::Vector{Int}
    field_size = length(apply_map)
    orders = zeros(Int, field_size)

    @inbounds for x in 0:(field_size - 1)
        orders[x + 1] != 0 && continue

        orbit = Int[]
        current = x

        while orders[current + 1] == 0 && !(current in orbit)
            push!(orbit, current)
            current = apply_map[current + 1]
        end

        orbit_length = length(orbit)
        for value in orbit
            orders[value + 1] = orbit_length
        end
    end

    return orders
end
