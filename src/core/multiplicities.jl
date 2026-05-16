using Nemo

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
    field = GF(2, n, "g")
    inputs = field_elements(field, n)
    function_values = function_values_to_field(lut, field, n)
    space_size = length(inputs)
    table = Matrix{Int}(undef, space_size, space_size)

    for (a_index, a) in pairs(inputs)
        for (b_index, b) in pairs(inputs)
            table[a_index, b_index] = walsh_coefficient(function_values, inputs, a, b)
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

    field = GF(2, n, "g")
    elements = field_elements(field, n)
    walsh_table = walsh_coefficient_table(lut, n)
    divisor = 2^(2 * n)
    multiplicities = Dict{Int, Int}()

    for (s_index, s) in pairs(elements)
        total = 0

        for a_index in eachindex(elements)
            for (b_index, b) in pairs(elements)
                total += trace_sign(b * s) * walsh_table[a_index, b_index]^k
            end
        end

        total % divisor == 0 || error("multiplicity sum is not divisible by $divisor")
        multiplicities[s_index - 1] = div(total, divisor)
    end

    return multiplicities
end
