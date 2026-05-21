abstract type APNCoefficient end

struct OneCoefficient <: APNCoefficient end

struct PowerCoefficient <: APNCoefficient
    base::Symbol
    exponent::Int
end

const ONE_COEFF = OneCoefficient()

struct APNTerm
    coefficient::APNCoefficient
    exponent::Int
end

struct APNTraceTerm
    n::Int
    step::Union{Nothing, Int}
    terms::Vector{APNTerm}
end

abstract type APNComponent end

struct APNFunction
    n::Union{Nothing, Int}
    id::Union{Nothing, String}
    terms::Vector{APNTerm}
    traces::Vector{APNTraceTerm}
    components::Vector{APNComponent}
end

struct APNFamilyMatch
    family::Symbol
    parameters::Dict{Symbol, Any}
    exact::Bool
    notes::Vector{String}
end

struct APNMonomial <: APNComponent
    coefficient_power::Int
    exponent::Int
end

struct APNTraceMonomial
    coefficient_power::Int
    exponent::Int
end

struct APNAbsoluteTrace <: APNComponent
    scale_power::Int
    terms::Vector{APNTraceMonomial}
end

struct APNRelativeTrace <: APNComponent
    scale_power::Int
    extension_degree::Int
    terms::Vector{APNTraceMonomial}
end

struct APNReference <: APNComponent
    id::String
end

struct Catalogue
    functions::Vector{APNFunction}
end

monomial_expr(exponent::Int) =
    APNTerm(ONE_COEFF, exponent)

monomial_expr(coefficient_power::Int, exponent::Int; base::Symbol = :p) =
    APNTerm(PowerCoefficient(base, coefficient_power), exponent)

function Tr(n::Int, terms::APNTerm...; step::Union{Nothing, Int} = nothing)
    return APNTraceTerm(n, step, collect(terms))
end

function APNFunction(n::Union{Nothing, Int}, id::Union{Nothing, AbstractString}, parts...)
    terms = APNTerm[]
    traces = APNTraceTerm[]
    components = APNComponent[]

    for part in parts
        if part isa APNTerm
            push!(terms, part)
        elseif part isa APNTraceTerm
            push!(traces, part)
        elseif part isa APNFunction
            append!(terms, part.terms)
            append!(traces, part.traces)
            append!(components, part.components)
        elseif part isa APNComponent
            push!(components, part)
        else
            throw(ArgumentError("unsupported APN function part: $part"))
        end
    end

    return APNFunction(n, id === nothing ? nothing : String(id), terms, traces, components)
end

APNFunction(n::Union{Nothing, Int}, parts...) = APNFunction(n, nothing, parts...)
APNFunction(id::AbstractString, parts...) = APNFunction(nothing, id, parts...)
APNFunction(parts...) = APNFunction(nothing, parts...)

Catalogue(functions::APNFunction...) = Catalogue(APNFunction[functions...])

trace_term(coefficient_power::Int, exponent::Int) =
    APNTraceMonomial(coefficient_power, exponent)

xpow(exponent::Int) =
    APNMonomial(0, exponent)

monomial(coefficient_power::Int, exponent::Int) =
    APNMonomial(coefficient_power, exponent)

reference(id::AbstractString) =
    APNReference(String(id))

function absolute_trace(terms::APNTraceMonomial...; scale::Int = 0)
    return APNAbsoluteTrace(scale, collect(terms))
end

function relative_trace(extension_degree::Int, terms::APNTraceMonomial...; scale::Int = 0)
    return APNRelativeTrace(scale, extension_degree, collect(terms))
end

Base.show(io::IO, ::OneCoefficient) = print(io, "1")

function Base.show(io::IO, coefficient::PowerCoefficient)
    print(io, coefficient.base)
    coefficient.exponent != 1 && print(io, coefficient.exponent)
end

function Base.show(io::IO, term::APNTerm)
    if term.coefficient isa OneCoefficient
        print(io, "x", term.exponent)
    else
        print(io, term.coefficient, "x", term.exponent)
    end
end

function Base.show(io::IO, trace::APNTraceTerm)
    print(io, "Tr", trace.n)
    trace.step !== nothing && print(io, "^", trace.step)
    print(io, "(")

    for (index, term) in enumerate(trace.terms)
        index > 1 && print(io, "+")
        print(io, term)
    end

    print(io, ")")
end

function Base.show(io::IO, monomial::APNMonomial)
    monomial.coefficient_power != 0 && print(io, "u", monomial.coefficient_power, "*")
    print(io, "x", monomial.exponent)
end

function Base.show(io::IO, monomial::APNTraceMonomial)
    monomial.coefficient_power != 0 && print(io, "u", monomial.coefficient_power, "*")
    print(io, "x", monomial.exponent)
end

function _show_trace_monomials(io::IO, terms::Vector{APNTraceMonomial})
    for (index, term) in enumerate(terms)
        index > 1 && print(io, "+")
        print(io, term)
    end
end

function Base.show(io::IO, trace::APNAbsoluteTrace)
    trace.scale_power != 0 && print(io, "u", trace.scale_power, "*")
    print(io, "tr(")
    _show_trace_monomials(io, trace.terms)
    print(io, ")")
end

function Base.show(io::IO, trace::APNRelativeTrace)
    trace.scale_power != 0 && print(io, "u", trace.scale_power, "*")
    print(io, "tr_", 2^trace.extension_degree, "/2(")
    _show_trace_monomials(io, trace.terms)
    print(io, ")")
end

function Base.show(io::IO, reference::APNReference)
    print(io, "(No. ", reference.id, ")")
end

function Base.show(io::IO, function_::APNFunction)
    parts = Any[function_.terms...]
    append!(parts, function_.traces)
    append!(parts, function_.components)
    isempty(parts) && return print(io, "0")

    for (index, part) in enumerate(parts)
        index > 1 && print(io, "+")
        print(io, part)
    end
end
