struct APNRepresentative
    id::Symbol
    family::Symbol
    n::Int
    function_::Union{Nothing, APNFunction}
    lut::Union{Nothing, Vector{Int}}
    lut_loader::Union{Nothing, Function}
    metadata::Dict{Symbol, Any}
end

struct APNEquivalenceResult
    representative::APNRepresentative
    equivalence::EAEquivalence
end

function APNRepresentative(id::Symbol,
                           family::Symbol,
                           n::Int,
                           function_::APNFunction;
                           metadata::AbstractDict = Dict{Symbol, Any}())
    return APNRepresentative(id,
                             family,
                             n,
                             apn_with_dimension(function_, n),
                             nothing,
                             nothing,
                             Dict{Symbol, Any}(metadata))
end

function APNRepresentative(id::Symbol,
                           family::Symbol,
                           n::Int,
                           lut::AbstractVector{<:Integer};
                           metadata::AbstractDict = Dict{Symbol, Any}())
    return APNRepresentative(id,
                             family,
                             n,
                             nothing,
                             lut_from_table(lut, n),
                             nothing,
                             Dict{Symbol, Any}(metadata))
end

function APNRepresentative(id::Symbol,
                           family::Symbol,
                           n::Int,
                           lut_loader::Function;
                           metadata::AbstractDict = Dict{Symbol, Any}())
    return APNRepresentative(id,
                             family,
                             n,
                             nothing,
                             nothing,
                             lut_loader,
                             Dict{Symbol, Any}(metadata))
end

function representative_lut(representative::APNRepresentative)::Vector{Int}
    representative.lut !== nothing && return copy(representative.lut)
    representative.function_ !== nothing && return apn_to_lut(representative.function_)
    representative.lut_loader !== nothing && return lut_from_table(representative.lut_loader(), representative.n)
    error("representative $(representative.id) has no function or LUT")
end

function identity_ea_equivalence(n::Int)::EAEquivalence
    identity = collect(0:(space_size(n) - 1))
    zero_map = zeros(Int, space_size(n))
    return EAEquivalence(identity, identity, zero_map)
end

function family_representatives(n::Int)::Vector{APNRepresentative}
    haskey(APN_REPRESENTATIVES_BY_DIMENSION, n) ||
        throw(ArgumentError("no representatives for dimension $n"))

    return copy(APN_REPRESENTATIVES_BY_DIMENSION[n])
end

function find_equivalences(n::Int,
                           lut::Union{AbstractVector{<:Integer}, AbstractDict{<:Integer, <:Integer}};
                           representatives::AbstractVector{<:APNRepresentative} = family_representatives(n),
                           k::Int = 4,
                           is_quadratic::Bool = true,
                           log_level::Union{Symbol, AbstractString} = :quiet,
                           max_external_maps_per_representative::Union{Nothing, Int} = nothing,
                           stop_at_first::Bool = true)::Vector{APNEquivalenceResult}
    if max_external_maps_per_representative !== nothing && max_external_maps_per_representative < 1
        throw(ArgumentError("max_external_maps_per_representative must be positive or nothing"))
    end

    target_lut = lut_from_table(lut, n)
    results = APNEquivalenceResult[]

    for representative in representatives
        representative.n == n ||
            throw(ArgumentError("representative $(representative.id) has dimension $(representative.n), expected $n"))

        if representative_lut(representative) == target_lut
            push!(results, APNEquivalenceResult(representative, identity_ea_equivalence(n)))
            stop_at_first && return results
        end
    end

    for representative in representatives
        equivalence = first_ea_equivalence(representative_lut(representative),
                                           target_lut,
                                           n,
                                           k = k,
                                           is_quadratic = is_quadratic,
                                           max_external_maps = max_external_maps_per_representative,
                                           log_level = log_level)
        equivalence === nothing && continue

        push!(results, APNEquivalenceResult(representative, equivalence))
        stop_at_first && return results
    end

    return results
end
