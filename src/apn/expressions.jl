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

struct APNFunction
    n::Union{Nothing, Int}
    terms::Vector{APNTerm}
    traces::Vector{APNTraceTerm}
end

struct APNFamilyMatch
    family::Symbol
    parameters::Dict{Symbol, Any}
    exact::Bool
    notes::Vector{String}
end

abstract type APNComponent end

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

struct APNDefinition
    n::Int
    table_id::String
    equation_id::String
    formula::String
    components::Vector{APNComponent}
end

x(exponent::Int) = APNTerm(ONE_COEFF, exponent)
pterm(power::Int, exponent::Int; base::Symbol = :p) =
    APNTerm(PowerCoefficient(base, power), exponent)
pterm(exponent::Int; base::Symbol = :p) = pterm(1, exponent; base = base)

function Tr(n::Int, terms::APNTerm...; step::Union{Nothing, Int} = nothing)
    return APNTraceTerm(n, step, collect(terms))
end

function APNFunction(n::Union{Nothing, Int}, parts...)
    terms = APNTerm[]
    traces = APNTraceTerm[]

    for part in parts
        if part isa APNTerm
            push!(terms, part)
        elseif part isa APNTraceTerm
            push!(traces, part)
        elseif part isa APNFunction
            append!(terms, part.terms)
            append!(traces, part.traces)
        else
            error("Unsupported APN function part: $part")
        end
    end

    return APNFunction(n, terms, traces)
end

APNFunction(parts...) = APNFunction(nothing, parts...)

trace_term(coefficient_power::Int, exponent::Int) =
    APNTraceMonomial(coefficient_power, exponent)
xpow(exponent::Int) = APNMonomial(0, exponent)
monomial(coefficient_power::Int, exponent::Int) =
    APNMonomial(coefficient_power, exponent)
reference(id::AbstractString) = APNReference(String(id))

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

function Base.show(io::IO, function_::APNFunction)
    parts = Any[function_.terms...]
    append!(parts, function_.traces)
    isempty(parts) && return print(io, "0")

    for (index, part) in enumerate(parts)
        index > 1 && print(io, "+")
        print(io, part)
    end
end
