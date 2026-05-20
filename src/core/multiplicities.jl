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
