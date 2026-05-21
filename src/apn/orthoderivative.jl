function _kernel_vector!(rows::Vector{Int}, pivot_cols::Vector{Int}, n::Int)::Int
    rank = 0

    @inbounds for col in 0:(n - 1)
        pivot = 0
        bit = 1 << col

        for row in (rank + 1):n
            if (rows[row] & bit) != 0
                pivot = row
                break
            end
        end

        pivot == 0 && continue

        rank += 1
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_cols[rank] = col

        pivot_row = rows[rank]
        for row in 1:n
            if row != rank && (rows[row] & bit) != 0
                rows[row] ⊻= pivot_row
            end
        end
    end

    rank < n || error("orthoderivative found a trivial kernel")

    free_col = 0
    @inbounds for col in 0:(n - 1)
        is_pivot = false

        for row in 1:rank
            if pivot_cols[row] == col
                is_pivot = true
                break
            end
        end

        if !is_pivot
            free_col = col
            break
        end
    end

    y = 1 << free_col

    @inbounds for row in rank:-1:1
        if isodd(count_ones(rows[row] & y))
            y |= 1 << pivot_cols[row]
        end
    end

    return y
end

@doc"""
    orthoderivative(F::Vector{Int}) -> Vector{Int}

Compute the orthoderivative map used by the MAGMA implementation for a
function represented by a lookup table over ``GF(2)^n``.

The lookup table is indexed by integers: the value of ``F(x)`` is stored at
`F[x + 1]`. The output coordinates are first computed as

    Trace(a^j * F(x)), j = 0, ..., n - 1

where `a` is the generator of ``GF(2)^n``. For each nonzero ``α``, the matrix
`J` has entries `J[j, i] = T_j(α ⊻ 2^i) + T_j(α) + T_j(2^i) + T_j(0)`.
The returned value is a generator of `Kernel(J)` encoded as an integer.
"""
function orthoderivative(F::Vector{Int})::Vector{Int}
    field_size = length(F)
    ispow2(field_size) || error("F length must be a power of two")

    n = trailing_zeros(field_size)
    check_lut_values(F, n, name = "F")

    field = GF(2, n, "g")
    generator = gen(field)
    function_values = function_values_to_field(F, field, n)
    trace_table = Matrix{Bool}(undef, n, field_size)

    @inbounds for component in 0:(n - 1)
        multiplier = generator^component

        for x in 0:(field_size - 1)
            trace_table[component + 1, x + 1] = isone(absolute_trace_bit(multiplier * function_values[x + 1]))
        end
    end

    pi = zeros(Int, field_size)
    rows = Vector{Int}(undef, n)
    pivot_cols = Vector{Int}(undef, n)

    @inbounds for alpha in 1:(field_size - 1)
        fill!(rows, 0)
        for component in 0:(n - 1)
            for i in 0:(n - 1)
                basis = 1 << i
                if trace_table[component + 1, (alpha ⊻ basis) + 1] ⊻
                   trace_table[component + 1, alpha + 1] ⊻
                   trace_table[component + 1, basis + 1] ⊻
                   trace_table[component + 1, 1]

                    rows[component + 1] |= (1 << i)
                end
            end
        end

        pi[alpha + 1] = _kernel_vector!(rows, pivot_cols, n)
    end

    return pi
end
