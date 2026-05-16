struct MultiplicityPartition
    blocks::Vector{Vector{Int}}
    block_index::Vector{Int}
    multiplicities::Vector{Int}
end

function lut_from_table(table::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    check_lut_values(table, n)
    return Int.(table)
end

function lut_from_table(table::AbstractDict{<:Integer, <:Integer}, n::Int)::Vector{Int}
    field_size = space_size(n)
    lut = fill(-1, field_size)

    for (x, y) in table
        check_space_value(x, n, name = "input values")
        check_space_value(y, n, name = "output values")
        lut[Int(x) + 1] == -1 || error("table contains a repeated input: $x")
        lut[Int(x) + 1] = Int(y)
    end

    all(value -> value != -1, lut) || error("table is missing at least one input")
    return lut
end

function partition_by_multiplicity(lut::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                   n::Int;
                                   k::Int = 4)::MultiplicityPartition
    normalized_lut = lut_from_table(lut, n)
    multiplicities_by_value = multiplicities_sigma(normalized_lut, n, k)
    values_by_multiplicity = Dict{Int, Vector{Int}}()

    for value in 0:(space_size(n) - 1)
        multiplicity = multiplicities_by_value[value]
        push!(get!(values_by_multiplicity, multiplicity, Int[]), value)
    end

    multiplicities = sort!(collect(keys(values_by_multiplicity)))
    blocks = Vector{Vector{Int}}(undef, length(multiplicities))
    block_index = zeros(Int, space_size(n))

    for (index, multiplicity) in pairs(multiplicities)
        block = sort!(values_by_multiplicity[multiplicity])
        blocks[index] = block

        for value in block
            block_index[value + 1] = index
        end
    end

    return MultiplicityPartition(blocks, block_index, multiplicities)
end

function gf2_rank(values::AbstractVector{<:Integer}, n::Int)::Int
    pivots = zeros(Int, n)
    rank_value = 0

    for raw_value in values
        value = Int(raw_value)

        while value != 0
            bit_index = 63 - leading_zeros(value)

            if pivots[bit_index + 1] == 0
                pivots[bit_index + 1] = value
                rank_value += 1
                break
            end

            value = xor(value, pivots[bit_index + 1])
        end
    end

    return rank_value
end

function gf2_is_independent(value::Int, pivots::Vector{Int})::Bool
    reduced = value

    while reduced != 0
        bit_index = 63 - leading_zeros(reduced)
        pivot = pivots[bit_index + 1]
        pivot == 0 && return true
        reduced = xor(reduced, pivot)
    end

    return false
end

function gf2_add_pivot!(pivots::Vector{Int}, value::Int)::Bool
    reduced = value

    while reduced != 0
        bit_index = 63 - leading_zeros(reduced)

        if pivots[bit_index + 1] == 0
            pivots[bit_index + 1] = reduced
            return true
        end

        reduced = xor(reduced, pivots[bit_index + 1])
    end

    return false
end

function gf2_basis(values::AbstractVector{<:Integer}, n::Int)::Vector{Int}
    pivots = zeros(Int, n)
    basis = Int[]

    for raw_value in values
        value = Int(raw_value)

        if gf2_add_pivot!(pivots, value)
            push!(basis, value)
            length(basis) == n && return basis
        end
    end

    return basis
end

function select_minimal_basis_union(partition::MultiplicityPartition, n::Int)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}}
    ordered_indices = sortperm(partition.blocks, by = length)
    block_count = length(ordered_indices)
    suffix_values = [Int[] for _ in 1:(block_count + 1)]

    for position in block_count:-1:1
        suffix_values[position] = vcat(partition.blocks[ordered_indices[position]], suffix_values[position + 1])
    end

    best_size = typemax(Int)
    best_indices = Int[]

    function search(position::Int, chosen::Vector{Int}, union_values::Vector{Int})::Nothing
        current_size = length(union_values)
        current_size >= best_size && return

        if gf2_rank(union_values, n) == n
            best_size = current_size
            empty!(best_indices)
            append!(best_indices, chosen)
            return
        end

        position > block_count && return
        gf2_rank(vcat(union_values, suffix_values[position]), n) == n || return

        block_index = ordered_indices[position]
        block = partition.blocks[block_index]

        push!(chosen, block_index)
        search(position + 1, chosen, vcat(union_values, block))
        pop!(chosen)

        search(position + 1, chosen, union_values)
        return
    end

    search(1, Int[], Int[])
    isempty(best_indices) && error("no union of multiplicity blocks contains a GF(2)-basis")

    union_values = sort!(reduce(vcat, (partition.blocks[index] for index in best_indices), init = Int[]))
    basis = gf2_basis(union_values, n)
    length(basis) == n || error("selected union does not contain a full basis")

    return basis, union_values, sort!(copy(best_indices))
end

function aligned_partitions(left::MultiplicityPartition, right::MultiplicityPartition)::Union{Nothing, MultiplicityPartition}
    length(left.blocks) == length(right.blocks) || return nothing

    right_by_multiplicity = Dict{Int, Vector{Int}}()
    for (index, multiplicity) in pairs(right.multiplicities)
        right_by_multiplicity[multiplicity] = right.blocks[index]
    end

    blocks = Vector{Vector{Int}}(undef, length(left.blocks))
    block_index = zeros(Int, length(right.block_index))

    for (index, multiplicity) in pairs(left.multiplicities)
        haskey(right_by_multiplicity, multiplicity) || return nothing
        block = right_by_multiplicity[multiplicity]
        length(block) == length(left.blocks[index]) || return nothing
        blocks[index] = block

        for value in block
            block_index[value + 1] = index
        end
    end

    return MultiplicityPartition(blocks, block_index, copy(left.multiplicities))
end

function backtrack_external_linear_maps(left::MultiplicityPartition,
                                        right::MultiplicityPartition,
                                        basis::Vector{Int},
                                        candidate_values::Vector{Int},
                                        n::Int)::Vector{Vector{Int}}
    check_length(basis, n, name = "basis")
    field_size = space_size(n)
    image = fill(-1, field_size)
    image[1] = 0
    domain_span = [0]
    image_pivots = zeros(Int, n)
    results = Vector{Vector{Int}}()

    _backtrack_external_linear_maps!(results,
                                     left,
                                     right,
                                     basis,
                                     candidate_values,
                                     n,
                                     image,
                                     domain_span,
                                     image_pivots,
                                     1)
    return results
end

const EA_EXTERNAL_LOG_LEVELS = Dict(:quiet => 0, :info => 1, :debug => 2)

function ea_external_log_level(level::Union{Symbol, AbstractString})::Int
    key = Symbol(level)
    haskey(EA_EXTERNAL_LOG_LEVELS, key) || error("log_level must be one of: quiet, info, debug")
    return EA_EXTERNAL_LOG_LEVELS[key]
end

function ea_external_log(level::Int, threshold::Int, message::String)::Nothing
    level <= threshold || return
    println(message)
    flush(stdout)
    return
end

function _backtrack_external_linear_maps!(results::Vector{Vector{Int}},
                                          left::MultiplicityPartition,
                                          right::MultiplicityPartition,
                                          basis::Vector{Int},
                                          candidate_values::Vector{Int},
                                          n::Int,
                                          image::Vector{Int},
                                          domain_span::Vector{Int},
                                          image_pivots::Vector{Int},
                                          start_level::Int)::Nothing
    function assign_basis_image(level::Int)::Nothing
        if level > n
            push!(results, copy(image))
            return
        end

        basis_value = basis[level]

        for candidate in candidate_values
            gf2_is_independent(candidate, image_pivots) || continue

            valid = true
            new_domain_values = Vector{Int}(undef, length(domain_span))
            new_image_values = Vector{Int}(undef, length(domain_span))

            @inbounds for index in eachindex(domain_span)
                x = domain_span[index]
                y = xor(x, basis_value)
                y_image = xor(image[x + 1], candidate)

                if right.block_index[y_image + 1] != left.block_index[y + 1]
                    valid = false
                    break
                end

                new_domain_values[index] = y
                new_image_values[index] = y_image
            end

            valid || continue

            old_span_length = length(domain_span)
            old_pivots = copy(image_pivots)
            gf2_add_pivot!(image_pivots, candidate)

            @inbounds for index in eachindex(new_domain_values)
                image[new_domain_values[index] + 1] = new_image_values[index]
                push!(domain_span, new_domain_values[index])
            end

            assign_basis_image(level + 1)

            resize!(domain_span, old_span_length)
            image_pivots .= old_pivots

            @inbounds for value in new_domain_values
                image[value + 1] = -1
            end
        end

        return
    end

    assign_basis_image(start_level)
    return
end

function backtrack_external_linear_maps_parallel(left::MultiplicityPartition,
                                                 right::MultiplicityPartition,
                                                 basis::Vector{Int},
                                                 candidate_values::Vector{Int},
                                                 n::Int;
                                                 log_level::Union{Symbol, AbstractString} = :quiet)::Vector{Vector{Int}}
    check_length(basis, n, name = "basis")
    field_size = space_size(n)
    isempty(candidate_values) && return Vector{Vector{Int}}()
    Threads.nthreads() > 1 || return backtrack_external_linear_maps(left, right, basis, candidate_values, n)

    verbosity = ea_external_log_level(log_level)
    basis_value = basis[1]
    branch_results = [Vector{Vector{Int}}() for _ in eachindex(candidate_values)]
    progress_lock = ReentrantLock()
    completed = Ref(0)
    progress_step = max(1, length(candidate_values) ÷ 10)

    ea_external_log(1, verbosity, "[info] external backtracking: $(length(candidate_values)) first-level branches on $(Threads.nthreads()) threads")

    function log_branch(candidate_index::Int, local_count::Int)::Nothing
        verbosity >= 1 || return

        lock(progress_lock)
        try
            completed[] += 1
            done = completed[]

            if verbosity >= 2
                println("[debug] external branch candidate_index=$candidate_index thread=$(Threads.threadid()) solutions=$local_count completed=$done/$(length(candidate_values))")
                flush(stdout)
            elseif done == 1 || done == length(candidate_values) || done % progress_step == 0
                println("[info] external backtracking progress completed=$done/$(length(candidate_values))")
                flush(stdout)
            end
        finally
            unlock(progress_lock)
        end

        return
    end

    Threads.@threads :static for candidate_index in eachindex(candidate_values)
        candidate = candidate_values[candidate_index]
        if candidate == 0
            log_branch(candidate_index, 0)
            continue
        end

        image = fill(-1, field_size)
        image[1] = 0
        domain_span = [0]
        image_pivots = zeros(Int, n)

        valid = true
        new_domain_values = Vector{Int}(undef, length(domain_span))
        new_image_values = Vector{Int}(undef, length(domain_span))

        @inbounds for index in eachindex(domain_span)
            x = domain_span[index]
            y = xor(x, basis_value)
            y_image = xor(image[x + 1], candidate)

            if right.block_index[y_image + 1] != left.block_index[y + 1]
                valid = false
                break
            end

            new_domain_values[index] = y
            new_image_values[index] = y_image
        end

        if !valid
            log_branch(candidate_index, 0)
            continue
        end

        gf2_add_pivot!(image_pivots, candidate)

        @inbounds for index in eachindex(new_domain_values)
            image[new_domain_values[index] + 1] = new_image_values[index]
            push!(domain_span, new_domain_values[index])
        end

        local_results = branch_results[candidate_index]
        _backtrack_external_linear_maps!(local_results,
                                         left,
                                         right,
                                         basis,
                                         candidate_values,
                                         n,
                                         image,
                                         domain_span,
                                         image_pivots,
                                         2)
        log_branch(candidate_index, length(local_results))
    end

    results = Vector{Vector{Int}}()
    for local_results in branch_results
        append!(results, local_results)
    end

    return results
end

function reconstruct_external_linear_maps(F::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          G::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}},
                                          n::Int;
                                          k::Int = 4,
                                          parallel::Bool = Threads.nthreads() > 1,
                                          log_level::Union{Symbol, AbstractString} = :quiet)::Vector{Vector{Int}}
    verbosity = ea_external_log_level(log_level)

    left = MultiplicityPartition(Vector{Vector{Int}}(), Int[], Int[])
    ea_external_log(1, verbosity, "[info] external reconstruction: partitioning F")
    left_time = @elapsed begin
        left = partition_by_multiplicity(F, n, k = k)
    end
    ea_external_log(1, verbosity, "[info] external reconstruction: partitioned F in $(round(left_time; digits = 6)) s")

    right = MultiplicityPartition(Vector{Vector{Int}}(), Int[], Int[])
    ea_external_log(1, verbosity, "[info] external reconstruction: partitioning G")
    right_time = @elapsed begin
        right = partition_by_multiplicity(G, n, k = k)
    end
    ea_external_log(1, verbosity, "[info] external reconstruction: partitioned G in $(round(right_time; digits = 6)) s")

    aligned_right = aligned_partitions(left, right)
    if aligned_right === nothing
        ea_external_log(1, verbosity, "[info] external reconstruction: multiplicity partitions are not aligned")
        return Vector{Vector{Int}}()
    end

    ea_external_log(1, verbosity, "[info] external reconstruction: selecting basis and candidate values")
    basis, _, chosen_indices = select_minimal_basis_union(left, n)
    candidate_values = sort!(reduce(vcat, (aligned_right.blocks[index] for index in chosen_indices), init = Int[]))
    ea_external_log(1, verbosity, "[info] external reconstruction: basis=$(basis), candidate_values=$(length(candidate_values)), parallel=$parallel")

    backtrack_time = 0.0
    results = Vector{Vector{Int}}()
    if parallel
        backtrack_time = @elapsed begin
            results = backtrack_external_linear_maps_parallel(left,
                                                              aligned_right,
                                                              basis,
                                                              candidate_values,
                                                              n,
                                                              log_level = log_level)
        end
    else
        ea_external_log(1, verbosity, "[info] external backtracking: running serial search")
        backtrack_time = @elapsed begin
            results = backtrack_external_linear_maps(left, aligned_right, basis, candidate_values, n)
        end
    end

    ea_external_log(1, verbosity, "[info] external reconstruction: backtracking finished in $(round(backtrack_time; digits = 6)) s with $(length(results)) candidates")
    return results
end
