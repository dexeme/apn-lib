using Nemo

mutable struct APNSearchContext
    n::Int
    space_size::Int
    ord_A::Vector{Int}
    ord_B::Vector{Int}
    apply_A::Vector{Int}
    apply_B::Vector{Int}
    visit_order::Vector{Int}
    used_outputs::BitVector
    solutions::Vector{Vector{Int}}
    max_solutions::Int
    on_solution::Function
    save_results::Bool
    class_index::Union{Nothing, Int}
    deadline::Union{Nothing, Float64}
    timed_out::Bool
    verify_apn_on_solution::Bool
end

function int_to_field_element(value::Integer, field, n::Int)
    0 <= value < 2^n || error("value must be between 0 and $(2^n - 1)")

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
    space_size = 2^n
    length(lut) == space_size || error("lut must have $space_size entries")

    field = GF(2, n, "g")
    polynomial_ring, _ = Nemo.polynomial_ring(field, "x")

    inputs = [int_to_field_element(value, field, n) for value in 0:(space_size - 1)]
    outputs = [int_to_field_element(lut[value + 1], field, n) for value in 0:(space_size - 1)]

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

function check_sbox_ddt_sizes(sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    space_size = length(sbox)
    check_square(ddt, name = "ddt")
    size(ddt) == (space_size, space_size) || error("ddt must be $space_size x $space_size")
    return true
end

function check_sbox_space_size(sbox::Vector{Int}, n::Int)::Bool
    expected_size = 2^n
    length(sbox) == expected_size || error("sbox must have $expected_size entries")
    return true
end

function even_hamming_weight_differences(space_size::Int)::Vector{Int}
    return [alpha for alpha in 1:(space_size - 1) if iseven(count_ones(alpha))]
end

function updateDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int}, delta::Int)::Bool
    delta in (-2, 2) || error("Unsupported DDT delta: $delta")

    c_value = sbox[c + 1]
    c_value != -1 || error("sbox[$c] must be assigned before updating the DDT")
    space_size = length(sbox)

    @inbounds for alpha in even_hamming_weight_differences(space_size)
        paired_x = xor(c, alpha)
        paired_value = sbox[paired_x + 1]
        paired_value == -1 && continue

        out_diff = xor(c_value, paired_value)
        ddt_value = ddt[alpha + 1, out_diff + 1] + delta
        ddt[alpha + 1, out_diff + 1] = ddt_value
        delta == 2 && ddt_value > 2 && return false
        delta == -2 && ddt_value == 2 && break
    end

    return true
end

function addDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    return updateDDTInformationUnchecked(c, sbox, ddt, 2)
end

function removeDDTInformationUnchecked(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    return updateDDTInformationUnchecked(c, sbox, ddt, -2)
end

function addDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return addDDTInformationUnchecked(c, sbox, ddt)
end

function removeDDTInformation(c::Int, sbox::Vector{Int}, ddt::Matrix{Int})::Bool
    check_sbox_ddt_sizes(sbox, ddt)
    return removeDDTInformationUnchecked(c, sbox, ddt)
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

function nextFreePosition(sbox::Vector{Int}, visit_order::Vector{Int})::Int
    @inbounds for value in visit_order
        sbox[value + 1] == -1 && return value
    end

    return -1
end

function standard_visit_order(n::Int)::Vector{Int}
    return collect(0:(2^n - 1))
end

function offset_visit_order(n::Int, offset::Int)::Vector{Int}
    space_size = 2^n
    return [mod(index + offset, space_size) for index in 0:(space_size - 1)]
end

function c_reference_visit_order(n::Int, class_index::Union{Nothing, Int} = nothing)::Vector{Int}
    if n == 7 && class_index == 24
        return offset_visit_order(n, 16)
    elseif n == 7 && class_index == 26
        return offset_visit_order(n, 8)
    end

    return standard_visit_order(n)
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
                          on_solution::Function = sbox -> nothing,
                          save_results::Bool = false,
                          class_index::Union{Nothing, Int} = nothing,
                          visit_order::Vector{Int} = c_reference_visit_order(n, class_index),
                          timeout_seconds::Union{Nothing, Real} = nothing,
                          verify_apn_on_solution::Bool = true)
    space_size = 2^n
    length(apply_A) == space_size || error("apply_A must have $space_size entries")
    length(apply_B) == space_size || error("apply_B must have $space_size entries")
    length(visit_order) == space_size || error("visit_order must have $space_size entries")
    !save_results || class_index !== nothing || error("class_index is required when save_results is true")
    deadline = timeout_seconds === nothing ? nothing : time() + Float64(timeout_seconds)

    return APNSearchContext(
        n,
        space_size,
        orbit_orders(apply_A),
        orbit_orders(apply_B),
        apply_A,
        apply_B,
        visit_order,
        falses(space_size),
        Vector{Vector{Int}}(),
        max_solutions,
        on_solution,
        save_results,
        class_index,
        deadline,
        false,
        verify_apn_on_solution,
    )
end

function APNSearchContext(n::Int;
                          ord_A::Function,
                          ord_B::Function,
                          apply_A::Function,
                          apply_B::Function,
                          max_solutions::Int = typemax(Int),
                          on_solution::Function = sbox -> nothing,
                          save_results::Bool = false,
                          class_index::Union{Nothing, Int} = nothing,
                          visit_order::Vector{Int} = c_reference_visit_order(n, class_index),
                          timeout_seconds::Union{Nothing, Real} = nothing,
                          verify_apn_on_solution::Bool = true)
    space_size = 2^n
    !save_results || class_index !== nothing || error("class_index is required when save_results is true")

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

    length(visit_order) == space_size || error("visit_order must have $space_size entries")
    deadline = timeout_seconds === nothing ? nothing : time() + Float64(timeout_seconds)

    return APNSearchContext(
        n,
        space_size,
        ord_A_values,
        ord_B_values,
        apply_A_values,
        apply_B_values,
        visit_order,
        falses(space_size),
        Vector{Vector{Int}}(),
        max_solutions,
        on_solution,
        save_results,
        class_index,
        deadline,
        false,
        verify_apn_on_solution,
    )
end

function seed_used_outputs!(context::APNSearchContext, sbox::Vector{Int})
    check_sbox_space_size(sbox, context.n)
    fill!(context.used_outputs, false)

    @inbounds for y in sbox
        y == -1 && continue
        context.used_outputs[y + 1] && error("Initial sbox is not injective")
        context.used_outputs[y + 1] = true
    end
end

function search_timed_out!(context::APNSearchContext)::Bool
    context.deadline === nothing && return false
    time() <= context.deadline && return false

    context.timed_out = true
    return true
end

function record_complete_solution!(context::APNSearchContext, sbox::Vector{Int})
    solution = copy(sbox)
    solution_is_apn = is_apn(solution)
    println("solution_is_apn = $solution_is_apn")

    if context.verify_apn_on_solution && !solution_is_apn
        error("Completed S-box is not APN")
    end

    push!(context.solutions, solution)
    context.on_solution(solution)

    if context.save_results
        save_search_result_constant(solution, context.n, context.class_index)
    end

    return context.solutions
end

function assign_orbit!(start_x::Int, start_y::Int, orbit_length::Int,
                       sbox::Vector{Int}, ddt::Matrix{Int},
                       context::APNSearchContext)::Tuple{Bool, Vector{Int}}
    current_x = start_x
    current_y = start_y
    assigned_inputs = Int[]

    for _ in 1:orbit_length
        if sbox[current_x + 1] != -1 || context.used_outputs[current_y + 1]
            return false, assigned_inputs
        end

        sbox[current_x + 1] = current_y
        context.used_outputs[current_y + 1] = true
        push!(assigned_inputs, current_x)

        if !addDDTInformationUnchecked(current_x, sbox, ddt)
            return false, assigned_inputs
        end

        current_x = context.apply_A[current_x + 1]
        current_y = context.apply_B[current_y + 1]
    end

    return true, assigned_inputs
end

function rollback_orbit!(assigned_inputs::Vector{Int}, sbox::Vector{Int},
                         ddt::Matrix{Int}, context::APNSearchContext)
    for assigned_x in Iterators.reverse(assigned_inputs)
        removeDDTInformationUnchecked(assigned_x, sbox, ddt)
        context.used_outputs[sbox[assigned_x + 1] + 1] = false
        sbox[assigned_x + 1] = -1
    end
end

function nextVal(depth::Int, sbox::Vector{Int}, ddt::Matrix{Int}, context::APNSearchContext)
    check_sbox_ddt_sizes(sbox, ddt)
    length(context.solutions) >= context.max_solutions && return context.solutions
    search_timed_out!(context) && return context.solutions

    if isComplete(sbox)
        return record_complete_solution!(context, sbox)
    end

    x = nextFreePosition(sbox, context.visit_order)
    x == -1 && return context.solutions
    ord_x = context.ord_A[x + 1]

    @inbounds for y in 0:(context.space_size - 1)
        length(context.solutions) >= context.max_solutions && break
        context.used_outputs[y + 1] && continue
        ord_x == context.ord_B[y + 1] || continue

        assignment_is_valid, assigned_inputs = assign_orbit!(x, y, ord_x, sbox, ddt, context)

        if assignment_is_valid
            nextVal(depth + 1, sbox, ddt, context)
        end

        rollback_orbit!(assigned_inputs, sbox, ddt, context)
    end

    return context.solutions
end

function nextVal(depth::Int, sbox::Vector{Int}, ddt::Matrix{Int};
                 ord_A::Function,
                 ord_B::Function,
                 apply_A::Function,
                 apply_B::Function,
                 max_solutions::Int = typemax(Int),
                 on_solution::Function = sbox -> nothing,
                 save_results::Bool = false,
                 class_index::Union{Nothing, Int} = nothing,
                 visit_order::Union{Nothing, Vector{Int}} = nothing,
                 timeout_seconds::Union{Nothing, Real} = nothing,
                 verify_apn_on_solution::Bool = true)
    n = trailing_zeros(length(sbox))
    2^n == length(sbox) || error("sbox length must be a power of 2")
    selected_visit_order = visit_order === nothing ? c_reference_visit_order(n, class_index) : visit_order

    context = APNSearchContext(
        n,
        ord_A = ord_A,
        ord_B = ord_B,
        apply_A = apply_A,
        apply_B = apply_B,
        max_solutions = max_solutions,
        on_solution = on_solution,
        save_results = save_results,
        class_index = class_index,
        visit_order = selected_visit_order,
        timeout_seconds = timeout_seconds,
        verify_apn_on_solution = verify_apn_on_solution,
    )
    seed_used_outputs!(context, sbox)

    return nextVal(depth, sbox, ddt, context)
end

function APNSearch(n::Int, A, B;
                   max_solutions::Int = typemax(Int),
                   on_solution::Function = sbox -> println(sbox),
                   sbox::Vector{Int} = fill(-1, 2^n),
                   ddt::Matrix{Int} = zeros(Int, 2^n, 2^n),
                   save_results::Bool = false,
                   class_index::Union{Nothing, Int} = nothing,
                   visit_order::Vector{Int} = c_reference_visit_order(n, class_index),
                   timeout_seconds::Union{Nothing, Real} = nothing,
                   seed_zero::Bool = true,
                   verify_apn_on_solution::Bool = true)
    check_sbox_space_size(sbox, n)
    check_sbox_ddt_sizes(sbox, ddt)

    if seed_zero
        sbox[1] in (-1, 0) || error("C reference search requires sbox[0] to be 0")
        if sbox[1] == -1
            sbox[1] = 0
        end
        addDDTInformation(0, sbox, ddt)
    end

    context = APNSearchContext(
        n,
        linear_map_lut(B, n),
        linear_map_lut(A, n),
        max_solutions = max_solutions,
        on_solution = on_solution,
        save_results = save_results,
        class_index = class_index,
        visit_order = visit_order,
        timeout_seconds = timeout_seconds,
        verify_apn_on_solution = verify_apn_on_solution,
    )
    seed_used_outputs!(context, sbox)

    return nextVal(0, sbox, ddt, context)
end

function APNSearchClasses(n::Int, class_indices = "all";
                          excluded_class_indices = Int[],
                          max_solutions::Int = 1,
                          on_solution::Function = (class_index, sbox) -> println("class $class_index: $sbox"),
                          save_results::Bool = true,
                          timeout_seconds::Union{Nothing, Real} = nothing,
                          seed_zero::Bool = true,
                          verify_apn_on_solution::Bool = true,
                          tuples_dir::String = default_tuples_dir())
    classes = normalize_precomputed_tuple_classes(
        n,
        class_indices,
        excluded_class_indices = excluded_class_indices,
        tuples_dir = tuples_dir,
    )
    results = Dict{Int, Vector{Vector{Int}}}()

    for class_index in classes
        A, B = precomputed_tuple_matrices(n, class_index, tuples_dir = tuples_dir)
        results[class_index] = APNSearch(
            n,
            A,
            B,
            max_solutions = max_solutions,
            on_solution = sbox -> on_solution(class_index, sbox),
            save_results = save_results,
            class_index = class_index,
            visit_order = c_reference_visit_order(n, class_index),
            timeout_seconds = timeout_seconds,
            seed_zero = seed_zero,
            verify_apn_on_solution = verify_apn_on_solution,
        )
    end

    return results
end
