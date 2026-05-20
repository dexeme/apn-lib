using Nemo

function check_binary_extension_field(field::FqField, n::Int)::Bool
    characteristic(field) == 2 || error("field must have characteristic 2")
    degree(field) == n || error("field degree must be n")
    return true
end

@doc"""
    absolute_trace_to_field(element::FqFieldElem) -> FqFieldElem

Return the absolute trace of `element` in its binary extension field.

### Input
- `element::FqFieldElem`: An element of `GF(2^n)`.

### Output
- `FqFieldElem`: The trace value in the prime subfield, represented in the
  same parent field as `element`.
"""
function absolute_trace_to_field(element::FqFieldElem)::FqFieldElem
    field = parent(element)
    n = degree(field)
    trace_value = zero(field)

    for i in 0:(n - 1)
        trace_value += element^(2^i)
    end

    return trace_value
end

function trace_to_field(element::FqFieldElem)::FqFieldElem
    return absolute_trace_to_field(element)
end

@doc"""
    absolute_trace_bit(element::FqFieldElem) -> Int

Return the absolute trace of `element` as an integer bit.

### Input
- `element::FqFieldElem`: An element of `GF(2^n)`.

### Output
- `Int`: `0` or `1`.
"""
function absolute_trace_bit(element::FqFieldElem)::Int
    trace_value = absolute_trace_to_field(element)
    iszero(trace_value) && return 0
    trace_value == one(parent(element)) || error("trace value must be in GF(2)")
    return 1
end

function trace_sign(element::FqFieldElem)::Int
    return absolute_trace_bit(element) == 0 ? 1 : -1
end

function field_elements(field::FqField, n::Int)::Vector{FqFieldElem}
    check_binary_extension_field(field, n)
    return [int_to_field_element(value, field, n) for value in 0:(2^n - 1)]
end

function function_values_to_field(lut::AbstractVector{<:Integer}, field::FqField, n::Int)::Vector{FqFieldElem}
    check_sbox_space_size(lut, n)
    return [int_to_field_element(value, field, n) for value in lut]
end

@doc"""
    walsh_coefficient(function_values, inputs, a, b) -> Int

Compute the Walsh coefficient `W_F(a, b)` from precomputed field inputs and
function values.

### Input
- `function_values::AbstractVector{<:FqFieldElem}`: Values `F(x)` over the input list.
- `inputs::AbstractVector{<:FqFieldElem}`: Elements `x` of `GF(2^n)`.
- `a::FqFieldElem`: Input mask.
- `b::FqFieldElem`: Output mask.

### Output
- `Int`: The Walsh coefficient.
"""
function walsh_coefficient(function_values::AbstractVector{<:FqFieldElem},
                           inputs::AbstractVector{<:FqFieldElem},
                           a::FqFieldElem,
                           b::FqFieldElem)::Int
    length(function_values) == length(inputs) || error("function_values and inputs must have the same length")
    parent(a) == parent(b) || error("a and b must belong to the same field")
    parent(a) == parent(first(inputs)) || error("inputs must belong to the same field as a and b")

    coefficient = 0
    @inbounds for index in eachindex(inputs)
        x = inputs[index]
        coefficient += trace_sign(b * function_values[index] + a * x)
    end

    return coefficient
end

@doc"""
    walsh_coefficient(lut, a, b, n) -> Int

Compute `W_F(a, b)` for an integer LUT representation of `F`.

### Input
- `lut::AbstractVector{<:Integer}`: S-box values indexed by integers `0:(2^n - 1)`.
- `a::FqFieldElem`: Input mask in `GF(2^n)`.
- `b::FqFieldElem`: Output mask in `GF(2^n)`.
- `n::Int`: Binary field extension degree.

### Output
- `Int`: The Walsh coefficient.
"""
function walsh_coefficient(lut::AbstractVector{<:Integer},
                           a::FqFieldElem,
                           b::FqFieldElem,
                           n::Int)::Int
    field = parent(a)
    parent(b) == field || error("a and b must belong to the same field")
    inputs = field_elements(field, n)
    function_values = function_values_to_field(lut, field, n)

    return walsh_coefficient(function_values, inputs, a, b)
end

@doc"""
    walsh_spectrum(lut, n) -> Vector{Int}

Compute the full Walsh spectrum of an `(n, n)` function represented by a LUT.

### Input
- `lut::AbstractVector{<:Integer}`: S-box values indexed by integers `0:(2^n - 1)`.
- `n::Int`: Binary field extension degree.

### Output
- `Vector{Int}`: All coefficients `W_F(a, b)` for `a, b in GF(2^n)`.
"""
function walsh_spectrum(lut::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    check_sbox_space_size(lut, n)
    field = GF(2, n, "g")
    inputs = field_elements(field, n)
    function_values = function_values_to_field(lut, field, n)
    spectrum = Int[]

    for a in inputs
        for b in inputs
            push!(spectrum, walsh_coefficient(function_values, inputs, a, b))
        end
    end

    return spectrum
end

function walsh_spectrum(function_::APNFunction, n::Int)::Vector{Int}
    return walsh_spectrum(apn_to_lut(function_, n), n)
end
