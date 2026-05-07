function is_apn(func_table::AbstractVector{<:Integer})::Bool
    field_size = length(func_table)

    for a in 1:(field_size - 1)
        diff_values = Set{Int}()

        for x in 0:(field_size - 1)
            result = Int(func_table[xor(x, a) + 1]) ⊻ Int(func_table[x + 1])
            result in diff_values && continue
            push!(diff_values, result)
        end

        length(diff_values) == div(field_size, 2) || return false
    end

    return true
end
