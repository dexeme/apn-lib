using Nemo

mutable struct APNSearchContext
    n::Int
    space_size::Int
    ord_A::Vector{Int}
    ord_B::Vector{Int}
    apply_A::Vector{Int}
    apply_B::Vector{Int}
    used_outputs::BitVector
    solutions::Vector{Vector{Int}}
    max_solutions::Int
    on_solution::Function
end

function check_sbox_ddt_sizes(sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    space_size = length(sbox)
    check_square(ddt, name = "ddt")
    size(ddt) == (space_size, space_size) || error("ddt must be $space_size x $space_size")
    return true
end

function updateDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int}, delta::Int)::Bool
    c_value = sbox[c + 1]
    c_value != -1 || error("sbox[$c] must be assigned before updating the DDT")
    space_size = length(sbox)
    reached_boundary_after_removal = false

    @inbounds for alpha in 1:(space_size - 1)
        isodd(count_ones(alpha)) && continue

        paired_x = xor(c, alpha)
        paired_value = sbox[paired_x + 1]
        paired_value == -1 && continue

        out_diff = xor(c_value, paired_value)
        ddt_value = ddt[alpha + 1, out_diff + 1] + delta
        ddt[alpha + 1, out_diff + 1] = ddt_value

        if delta > 0
            ddt_value > 2 && return false
        elseif ddt_value == 2
            reached_boundary_after_removal = true
        end
    end

    return delta > 0 || !reached_boundary_after_removal
end

function addDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return updateDDTInformationUnchecked(c, sbox, ddt, 2)
end

function removeDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return updateDDTInformationUnchecked(c, sbox, ddt, -2)
end

function isComplete(sbox::Vector{Int})::Bool
    @inbounds for value in sbox
        value == -1 && return false
    end

    return true
end

function nextFreePosition(sbox::Vector{Int})::Int
    @inbounds for index in eachindex(sbox)
        sbox[index] == -1 && return index - 1
    end

    return -1
end

function int_matrix_to_lut(M::AbstractMatrix, n::Int)::Vector{Int}
    check_square(M)
    size(M) == (n, n) || error("Matrix must be $n x $n")
    space_size = 2^n
    table = Vector{Int}(undef, space_size)

    @inbounds for x in 0:(space_size - 1)
        output = 0

        for row in 1:n
            output_bit = false

            for col in 1:n
                if isodd((x >> (col - 1)) & 1) && isodd(Int(M[row, col]))
                    output_bit = !output_bit
                end
            end

            if output_bit
                output |= 1 << (row - 1)
            end
        end

        table[x + 1] = output
    end

    return table
end

function linear_map_lut(M::FqMatrix, n::Int)::Vector{Int}
    _ensure_gf2_matrix(M, n)
    return matrix_to_sbox(M)
end

function linear_map_lut(M::AbstractMatrix, n::Int)::Vector{Int}
    return int_matrix_to_lut(M, n)
end

function orbit_orders(apply_map::Vector{Int})::Vector{Int}
    space_size = length(apply_map)
    orders = zeros(Int, space_size)

    @inbounds for x in 0:(space_size - 1)
        orders[x + 1] != 0 && continue

        orbit = Int[]
        current = x

        while orders[current + 1] == 0 && !(current in orbit)
            push!(orbit, current)
            current = apply_map[current + 1]
        end

        orbit_length = length(orbit)
        for value in orbit
            orders[value + 1] = orbit_length
        end
    end

    return orders
end

function APNSearchContext(n::Int, apply_A::Vector{Int}, apply_B::Vector{Int};
                          max_solutions::Int = typemax(Int),
                          on_solution::Function = sbox -> nothing)
    space_size = 2^n
    length(apply_A) == space_size || error("apply_A must have $space_size entries")
    length(apply_B) == space_size || error("apply_B must have $space_size entries")

    return APNSearchContext(
        n,
        space_size,
        orbit_orders(apply_A),
        orbit_orders(apply_B),
        apply_A,
        apply_B,
        falses(space_size),
        Vector{Vector{Int}}(),
        max_solutions,
        on_solution,
    )
end

function APNSearchContext(n::Int;
                          ord_A::Function,
                          ord_B::Function,
                          apply_A::Function,
                          apply_B::Function,
                          max_solutions::Int = typemax(Int),
                          on_solution::Function = sbox -> nothing)
    space_size = 2^n
    apply_A_values = Vector{Int}(undef, space_size)
    apply_B_values = Vector{Int}(undef, space_size)
    ord_A_values = Vector{Int}(undef, space_size)
    ord_B_values = Vector{Int}(undef, space_size)

    @inbounds for x in 0:(space_size - 1)
        apply_A_values[x + 1] = apply_A(x)
        apply_B_values[x + 1] = apply_B(x)
        ord_A_values[x + 1] = ord_A(x)
        ord_B_values[x + 1] = ord_B(x)
    end

    return APNSearchContext(
        n,
        space_size,
        ord_A_values,
        ord_B_values,
        apply_A_values,
        apply_B_values,
        falses(space_size),
        Vector{Vector{Int}}(),
        max_solutions,
        on_solution,
    )
end

function seed_used_outputs!(context::APNSearchContext, sbox::Vector{Int})
    length(sbox) == context.space_size || error("sbox must have $(context.space_size) entries")
    fill!(context.used_outputs, false)

    @inbounds for y in sbox
        y == -1 && continue
        context.used_outputs[y + 1] && error("Initial sbox is not injective")
        context.used_outputs[y + 1] = true
    end
end

function nextVal(depth::Int, sbox::Vector{Int}, ddt::Matrix{Int}, context::APNSearchContext)
    check_sbox_ddt_sizes(sbox, ddt)
    length(context.solutions) >= context.max_solutions && return context.solutions

    if isComplete(sbox)
        solution = copy(sbox)
        push!(context.solutions, solution)
        context.on_solution(solution)
        return context.solutions
    end

    x = nextFreePosition(sbox)
    x == -1 && return context.solutions
    ord_x = context.ord_A[x + 1]

    @inbounds for y in 0:(context.space_size - 1)
        length(context.solutions) >= context.max_solutions && break
        context.used_outputs[y + 1] && continue
        ord_x == context.ord_B[y + 1] || continue

        current_x = x
        current_y = y
        assigned_inputs = Int[]
        assignment_is_valid = true

        for _ in 1:ord_x
            if sbox[current_x + 1] != -1 || context.used_outputs[current_y + 1]
                assignment_is_valid = false
                break
            end

            sbox[current_x + 1] = current_y
            context.used_outputs[current_y + 1] = true
            push!(assigned_inputs, current_x)

            if !updateDDTInformationUnchecked(current_x, sbox, ddt, 2)
                assignment_is_valid = false
                break
            end

            current_x = context.apply_A[current_x + 1]
            current_y = context.apply_B[current_y + 1]
        end

        if assignment_is_valid
            nextVal(depth + 1, sbox, ddt, context)
        end

        for assigned_x in Iterators.reverse(assigned_inputs)
            updateDDTInformationUnchecked(assigned_x, sbox, ddt, -2)
            context.used_outputs[sbox[assigned_x + 1] + 1] = false
            sbox[assigned_x + 1] = -1
        end
    end

    return context.solutions
end

function nextVal(depth::Int, sbox::Vector{Int}, ddt::Matrix{Int};
                 ord_A::Function,
                 ord_B::Function,
                 apply_A::Function,
                 apply_B::Function,
                 max_solutions::Int = typemax(Int),
                 on_solution::Function = sbox -> nothing)
    n = trailing_zeros(length(sbox))
    2^n == length(sbox) || error("sbox length must be a power of 2")

    context = APNSearchContext(
        n,
        ord_A = ord_A,
        ord_B = ord_B,
        apply_A = apply_A,
        apply_B = apply_B,
        max_solutions = max_solutions,
        on_solution = on_solution,
    )
    seed_used_outputs!(context, sbox)

    return nextVal(depth, sbox, ddt, context)
end

function APNSearch(n::Int, A, B;
                   max_solutions::Int = typemax(Int),
                   on_solution::Function = sbox -> println(sbox),
                   sbox::Vector{Int} = fill(-1, 2^n),
                   ddt::Matrix{Int} = zeros(Int, 2^n, 2^n))
    space_size = 2^n
    length(sbox) == space_size || error("sbox must have $space_size entries")
    size(ddt) == (space_size, space_size) || error("ddt must be $space_size x $space_size")

    context = APNSearchContext(
        n,
        linear_map_lut(A, n),
        linear_map_lut(B, n),
        max_solutions = max_solutions,
        on_solution = on_solution,
    )
    seed_used_outputs!(context, sbox)

    return nextVal(0, sbox, ddt, context)
end

function APNsearch(args...; kwargs...)
    return APNSearch(args...; kwargs...)
end
