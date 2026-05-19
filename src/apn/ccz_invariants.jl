@inline function _encoded_graph_point(x::Int, y::Int, n::Int)::Int
    return x | (y << n)
end

function _graph_points(lut::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    check_lut_values(lut, n)

    field_size = space_size(n)
    points = Vector{Int}(undef, field_size)

    @inbounds for x in 0:(field_size - 1)
        y = Int(lut[x + 1])
        points[x + 1] = _encoded_graph_point(x, y, n)
    end

    return points
end

function _delta_points(lut::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    check_lut_values(lut, n)

    field_size = space_size(n)
    points = Set{Int}()

    @inbounds for a in 1:(field_size - 1)
        for x in 0:(field_size - 1)
            y = Int(lut[(x ⊻ a) + 1]) ⊻ Int(lut[x + 1])
            point = _encoded_graph_point(a, y, n)
            push!(points, point)
        end
    end

    return sort!(collect(points))
end

@inline function _set_bit!(words::Vector{UInt64}, bit_index::Int)::Vector{UInt64}
    word_index = div(bit_index, 64) + 1
    offset = bit_index % 64

    @inbounds words[word_index] |= UInt64(1) << offset

    return words
end

function _translated_indicator(points::AbstractVector{<:Integer},
                               translation::Int,
                               word_count::Int)::Vector{UInt64}
    words = zeros(UInt64, word_count)

    @inbounds for point in points
        translated = Int(point) ⊻ translation
        _set_bit!(words, translated)
    end

    return words
end

function _leading_bit_index(words::Vector{UInt64})::Int
    @inbounds for word_index in length(words):-1:1
        word = words[word_index]

        if word != 0
            return (word_index - 1) * 64 + (63 - leading_zeros(word))
        end
    end

    return -1
end

function _xor_words!(target::Vector{UInt64}, pivot::Vector{UInt64})::Vector{UInt64}
    @inbounds for i in eachindex(target)
        target[i] ⊻= pivot[i]
    end

    return target
end

function _development_rank(points::AbstractVector{<:Integer}, ambient_dim::Int)::Int
    ambient_size = space_size(ambient_dim)
    word_count = cld(ambient_size, 64)

    pivots = Dict{Int, Vector{UInt64}}()
    rank_value = 0

    for translation in 0:(ambient_size - 1)
        row = _translated_indicator(points, translation, word_count)

        while true
            pivot_bit = _leading_bit_index(row)

            if pivot_bit == -1
                break
            end

            if haskey(pivots, pivot_bit)
                _xor_words!(row, pivots[pivot_bit])
            else
                pivots[pivot_bit] = row
                rank_value += 1
                break
            end
        end
    end

    return rank_value
end

@doc"""
    gamma_rank(lut, n) -> Int

Compute the Γ-rank of a function `F : GF(2)^n -> GF(2)^n` represented by a LUT.

The Γ-rank is the rank over `GF(2)` of the development matrix of the graph

    Γ_F = {(x, F(x)) : x ∈ GF(2)^n}.

The graph point `(x, F(x))` is encoded as the integer

    x | (F(x) << n).
"""
function gamma_rank(lut::AbstractVector{<:Integer}, n::Int)::Int
    graph = _graph_points(lut, n)
    return _development_rank(graph, 2 * n)
end

@doc"""
    delta_rank(lut, n) -> Int

Compute the Δ-rank of an APN function `F : GF(2)^n -> GF(2)^n` represented by a LUT.

The Δ-rank is computed from the development of the differential set

    Δ_F = {(a, F(x + a) + F(x)) : a ≠ 0, x ∈ GF(2)^n}.

Since the field has characteristic two, both additions are implemented as XOR.
"""
function delta_rank(lut::AbstractVector{<:Integer}, n::Int)::Int
    delta = _delta_points(lut, n)
    return _development_rank(delta, 2 * n)
end