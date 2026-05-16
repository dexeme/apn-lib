using Nemo

struct ANFCoordinate
    n::Int
    coefficients::BitVector
end

struct ANFVector
    n::Int
    coordinates::Vector{ANFCoordinate}
end

function field_element_lookup(field::FqField, n::Int)::Dict{FqFieldElem, Int}
    elements = field_elements(field, n)
    return Dict(element => index - 1 for (index, element) in pairs(elements))
end

function field_element_to_int(element::FqFieldElem, n::Int)::Int
    field = parent(element)
    lookup = field_element_lookup(field, n)
    haskey(lookup, element) || error("element does not belong to GF(2^$n)")
    return lookup[element]
end

@doc"""
    univariate_to_lut(polynomial::FqPolyRingElem, n::Int) -> Vector{Int}

Evaluate a univariate polynomial over `GF(2^n)` on every field element.

### Input
- `polynomial::FqPolyRingElem`: Polynomial `F(x)` over `GF(2^n)`.
- `n::Int`: Binary field extension degree.

### Output
- `Vector{Int}`: Lookup table where index `x + 1` stores `F(x)` as an integer.
"""
function univariate_to_lut(polynomial::FqPolyRingElem, n::Int)::Vector{Int}
    field = base_ring(parent(polynomial))
    check_binary_extension_field(field, n)
    lookup = field_element_lookup(field, n)
    lut = Vector{Int}(undef, 2^n)

    for x_int in 0:(2^n - 1)
        x = int_to_field_element(x_int, field, n)
        y = polynomial(x)
        haskey(lookup, y) || error("polynomial output does not belong to GF(2^$n)")
        lut[x_int + 1] = lookup[y]
    end

    return lut
end

@doc"""
    lut_to_univariate(lut::AbstractVector{<:Integer}, n::Int) -> Tuple{FqPolyRingElem, FqField}

Interpolate the unique reduced univariate polynomial over `GF(2^n)` matching a
lookup table.

### Input
- `lut::AbstractVector{<:Integer}`: Function values indexed by integers `0:(2^n - 1)`.
- `n::Int`: Binary field extension degree.

### Output
- `Tuple{FqPolyRingElem, FqField}`: The interpolated polynomial and its field.
"""
function lut_to_univariate(lut::AbstractVector{<:Integer}, n::Int)
    check_lut_values(lut, n)
    return interpolate_sbox_polynomial(Int.(lut), n)
end

function truth_table_to_anf_coefficients(truth_table::AbstractVector{Bool}, n::Int)::BitVector
    check_space_length(truth_table, n, name = "truth table")
    coefficients = BitVector(truth_table)

    # Fast Möbius transform over GF(2): after this loop, each mask stores the
    # coefficient of the monomial selected by that mask.
    for bit in 0:(n - 1)
        step = 1 << bit
        for mask in 0:(2^n - 1)
            if (mask & step) != 0
                coefficients[mask + 1] = xor(coefficients[mask + 1], coefficients[(mask ⊻ step) + 1])
            end
        end
    end

    return coefficients
end

function anf_coefficients_to_truth_table(coefficients::BitVector, n::Int)::BitVector
    check_space_length(coefficients, n, name = "ANF coefficient vector")
    truth_table = copy(coefficients)

    # The Möbius transform is self-inverse over GF(2).
    for bit in 0:(n - 1)
        step = 1 << bit
        for mask in 0:(2^n - 1)
            if (mask & step) != 0
                truth_table[mask + 1] = xor(truth_table[mask + 1], truth_table[(mask ⊻ step) + 1])
            end
        end
    end

    return truth_table
end

@doc"""
    lut_to_anf(lut::AbstractVector{<:Integer}, n::Int) -> ANFVector

Convert a vectorial lookup table to Algebraic Normal Form using the fast
Möbius transform, one coordinate bit at a time.

### Input
- `lut::AbstractVector{<:Integer}`: Function values indexed by integers `0:(2^n - 1)`.
- `n::Int`: Number of input and output bits.

### Output
- `ANFVector`: Vectorial ANF with one `ANFCoordinate` for each output bit.
"""
function lut_to_anf(lut::AbstractVector{<:Integer}, n::Int)::ANFVector
    check_lut_values(lut, n)
    coordinates = ANFCoordinate[]

    for output_bit in 0:(n - 1)
        truth_table = [isodd((value >> output_bit) & 1) for value in lut]
        push!(coordinates, ANFCoordinate(n, truth_table_to_anf_coefficients(truth_table, n)))
    end

    return ANFVector(n, coordinates)
end

@doc"""
    anf_to_lut(anf::ANFVector) -> Vector{Int}

Evaluate a vectorial ANF on every input vector.

### Input
- `anf::ANFVector`: Vectorial Algebraic Normal Form.

### Output
- `Vector{Int}`: Lookup table where index `x + 1` stores `F(x)` as an integer.
"""
function anf_to_lut(anf::ANFVector)::Vector{Int}
    n = anf.n
    length(anf.coordinates) == n || error("ANF must have n coordinate functions")
    lut = zeros(Int, space_size(n))

    for (output_index, coordinate) in pairs(anf.coordinates)
        coordinate.n == n || error("all ANF coordinates must have degree n")
        truth_table = anf_coefficients_to_truth_table(coordinate.coefficients, n)

        for x_int in 0:(space_size(n) - 1)
            if truth_table[x_int + 1]
                lut[x_int + 1] |= 1 << (output_index - 1)
            end
        end
    end

    return lut
end

function lut_to_graph(lut::AbstractVector{<:Integer}, n::Int)::Vector{Tuple{Int, Int}}
    check_lut_values(lut, n)
    return [(x, Int(lut[x + 1])) for x in 0:(space_size(n) - 1)]
end

function graph_to_lut(graph::AbstractVector{<:Tuple{<:Integer, <:Integer}}, n::Int)::Vector{Int}
    field_size = space_size(n)
    check_space_length(graph, n, name = "graph", unit = "pairs")
    lut = fill(-1, field_size)

    for (x, y) in graph
        check_space_value(x, n, name = "graph input values")
        check_space_value(y, n, name = "graph output values")
        lut[Int(x) + 1] == -1 || error("graph contains a repeated input: $x")
        lut[Int(x) + 1] = Int(y)
    end

    all(value -> value != -1, lut) || error("graph is missing at least one input")
    return lut
end

univariate_to_graph(polynomial::FqPolyRingElem, n::Int)::Vector{Tuple{Int, Int}} =
    lut_to_graph(univariate_to_lut(polynomial, n), n)

anf_to_graph(anf::ANFVector)::Vector{Tuple{Int, Int}} =
    lut_to_graph(anf_to_lut(anf), anf.n)

graph_to_univariate(graph::AbstractVector{<:Tuple{<:Integer, <:Integer}}, n::Int) =
    lut_to_univariate(graph_to_lut(graph, n), n)

graph_to_anf(graph::AbstractVector{<:Tuple{<:Integer, <:Integer}}, n::Int)::ANFVector =
    lut_to_anf(graph_to_lut(graph, n), n)

univariate_to_anf(polynomial::FqPolyRingElem, n::Int)::ANFVector =
    lut_to_anf(univariate_to_lut(polynomial, n), n)

anf_to_univariate(anf::ANFVector) =
    lut_to_univariate(anf_to_lut(anf), anf.n)
