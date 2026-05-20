using Nemo

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
    multiplicities_sigma(lut, n, k) -> Dict{Int, Int}

Compute `M_F^k(0, s)` for every `s in GF(2^n)` using the Walsh-coefficient
formula.

### Input
- `lut::AbstractVector{<:Integer}`: S-box values indexed by integers `0:(2^n - 1)`.
- `n::Int`: Binary field extension degree.
- `k::Int`: Power used in the multiplicity formula.

### Output
- `Dict{Int, Int}`: Multiplicity by integer representation of `s`.
"""
function multiplicities_sigma(lut::AbstractVector{<:Integer}, n::Int, k::Int)::Dict{Int, Int}
    k >= 0 || error("k must be non-negative")
    check_sbox_space_size(lut, n)

    field_size = space_size(n)
    normalized_lut = Int.(lut)
    dual_indices = trace_dual_indices(n)
    transformed_sums = zeros(Int, field_size)
    divisor = field_size^2
    multiplicities = Dict{Int, Int}()

    for b in 0:(field_size - 1)
        output_mask = dual_indices[b + 1]
        column = Vector{Int}(undef, field_size)

        @inbounds for x in 0:(field_size - 1)
            column[x + 1] = parity_sign(output_mask & normalized_lut[x + 1])
        end

        fwht!(column)
        sum_over_a = 0

        @inbounds for coefficient in column
            sum_over_a += coefficient^k
        end

        transformed_sums[output_mask + 1] = sum_over_a
    end

    fwht!(transformed_sums)

    for s in 0:(field_size - 1)
        total = transformed_sums[s + 1]
        total % divisor == 0 || error("multiplicity sum is not divisible by $divisor")
        multiplicities[s] = div(total, divisor)
    end

    return multiplicities
end

function multiplicities_sigma(function_::APNFunction, n::Int, k::Int)::Dict{Int, Int}
    return multiplicities_sigma(apn_to_lut(function_, n), n, k)
end
