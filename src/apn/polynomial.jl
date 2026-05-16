using Nemo

function int_to_field_element(value::Integer, field, n::Int)
    check_space_value(value, n)

    coefficients = [isodd((value >> bit_index) & 1) ? 1 : 0 for bit_index in 0:(n - 1)]
    return field(coefficients)
end

function field_power_lookup(generator, n::Int)
    field = parent(generator)
    powers = Dict{typeof(generator), Int}()
    current = one(field)

    for exponent in 0:(2^n - 2)
        powers[current] = exponent
        current = current * generator
    end

    return powers
end

function interpolate_sbox_polynomial(lut::Vector{Int}, n::Int)
    check_lut_values(lut, n)
    field_size = space_size(n)

    field = GF(2, n, "g")
    polynomial_ring, _ = Nemo.polynomial_ring(field, "x")

    inputs = [int_to_field_element(value, field, n) for value in 0:(field_size - 1)]
    outputs = [int_to_field_element(lut[value + 1], field, n) for value in 0:(field_size - 1)]

    return interpolate(polynomial_ring, inputs, outputs), field
end

function format_sbox_polynomial(lut::Vector{Int}, n::Int)::String
    polynomial, field = interpolate_sbox_polynomial(lut, n)
    iszero(polynomial) && return "0"

    generator = gen(field)
    powers = field_power_lookup(generator, n)
    terms = String[]

    for exponent in 0:degree(polynomial)
        coefficient = coeff(polynomial, exponent)
        iszero(coefficient) && continue

        if coefficient == one(field)
            push!(terms, "x^$exponent")
        else
            generator_exponent = powers[coefficient]
            push!(terms, "g^$(generator_exponent)x^$exponent")
        end
    end

    return isempty(terms) ? "0" : join(terms, " + ")
end
