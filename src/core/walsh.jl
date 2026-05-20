function fwht!(values::Vector{Int})::Vector{Int}
    length(values) > 0 || error("FWHT input must be non-empty")
    ispow2(length(values)) || error("FWHT input length must be a power of two")

    half_block = 1
    value_count = length(values)

    while half_block < value_count
        block = 2 * half_block

        for start in 1:block:value_count
            @inbounds for index in start:(start + half_block - 1)
                x = values[index]
                y = values[index + half_block]
                values[index] = x + y
                values[index + half_block] = x - y
            end
        end

        half_block = block
    end

    return values
end

@inline function parity_sign(mask::Int)::Int
    return isodd(count_ones(mask)) ? -1 : 1
end

function trace_dual_indices(n::Int)::Vector{Int}
    field_size = space_size(n)
    field = GF(2, n, "g")
    elements = field_elements(field, n)
    dual_indices = Vector{Int}(undef, field_size)

    @inbounds for c in 0:(field_size - 1)
        c_element = elements[c + 1]
        dual_index = 0

        for bit in 0:(n - 1)
            basis_element = elements[(1 << bit) + 1]
            if absolute_trace_bit(c_element * basis_element) == 1
                dual_index |= 1 << bit
            end
        end

        dual_indices[c + 1] = dual_index
    end

    return dual_indices
end

@doc"""
    walsh_coefficient_table(lut, n) -> Matrix{Int}

Compute all Walsh coefficients of an `(n, n)` function represented by a LUT.

### Input
- `lut::AbstractVector{<:Integer}`: S-box values indexed by integers `0:(2^n - 1)`.
- `n::Int`: Binary field extension degree.

### Output
- `Matrix{Int}`: Table whose `(a + 1, b + 1)` entry stores `W_F(a, b)`.
"""
function walsh_coefficient_table(lut::AbstractVector{<:Integer}, n::Int)::Matrix{Int}
    check_sbox_space_size(lut, n)
    field_size = space_size(n)
    normalized_lut = Int.(lut)
    dual_indices = trace_dual_indices(n)
    table = Matrix{Int}(undef, field_size, field_size)

    for b in 0:(field_size - 1)
        output_mask = dual_indices[b + 1]
        column = Vector{Int}(undef, field_size)

        @inbounds for x in 0:(field_size - 1)
            column[x + 1] = parity_sign(output_mask & normalized_lut[x + 1])
        end

        fwht!(column)

        @inbounds for a in 0:(field_size - 1)
            table[a + 1, b + 1] = column[dual_indices[a + 1] + 1]
        end
    end

    return table
end

@doc"""
    walsh_spectrum(lut::AbstractVector{<:Integer}, n::Int) -> Vector{Int}

Return the pure Walsh spectrum of a function as a vector of length `2^(2n)`.
This function uses the coefficient table computed quickly via FWHT.
"""
function walsh_spectrum(lut::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    table = walsh_coefficient_table(lut, n)
    return vec(table)
end

function walsh_spectrum(function_::APNFunction, n::Int)::Vector{Int}
    return walsh_spectrum(apn_to_lut(function_, n), n)
end

@doc"""
    extended_walsh_spectrum(lut::AbstractVector{<:Integer}, n::Int) -> Vector{Int}

Return the extended Walsh spectrum, defined as the sorted multiset of absolute
Walsh coefficient values. This is useful for testing linear equivalence of
orthoderivatives.
"""
function extended_walsh_spectrum(lut::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    raw_spectrum = walsh_spectrum(lut, n)
    return sort(abs.(raw_spectrum))
end

function extended_walsh_spectrum(function_::APNFunction, n::Int)::Vector{Int}
    return extended_walsh_spectrum(apn_to_lut(function_, n), n)
end
