using Nemo

# Families C1-C11 from:
# https://boolean.w.uib.no/files/2018/07/quadratic_APN_poly.pdf
#
# The code below is intentionally symbolic.  It is meant to answer questions of
# the form "does this written APN polynomial have the shape of family C_i?"
# before one starts evaluating the polynomial over a concrete finite field.

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

function coefficient_key(coefficient::APNCoefficient)
    coefficient isa OneCoefficient && return (:one, 0)
    coefficient isa PowerCoefficient && return (coefficient.base, coefficient.exponent)
    error("Unknown coefficient type: $(typeof(coefficient))")
end

function term_key(term::APNTerm, n::Union{Nothing, Int})
    modulus = n === nothing ? nothing : 2^n - 1
    exponent = modulus === nothing ? term.exponent : normalize_exponent(term.exponent, modulus)
    return (coefficient_key(term.coefficient), exponent)
end

normalize_exponent(exponent::Int, modulus::Int) = mod(exponent - 1, modulus) + 1

function infer_n(function_::APNFunction)
    function_.n !== nothing && return function_.n

    trace_degrees = unique(trace.n for trace in function_.traces)
    length(trace_degrees) == 1 && return only(trace_degrees)
    return nothing
end

function expanded_terms(function_::APNFunction)
    n = infer_n(function_)
    terms = copy(function_.terms)
    n === nothing && return terms

    modulus = 2^n - 1
    for trace in function_.traces
        trace.n == n || error("Trace degree Tr$(trace.n) does not match function degree n=$n")
        step = trace.step === nothing ? 1 : trace.step
        iterations = div(n, step)

        for j in 0:(iterations - 1)
            multiplier = 2^(step * j)
            for term in trace.terms
                push!(terms, APNTerm(term.coefficient, normalize_exponent(term.exponent * multiplier, modulus)))
            end
        end
    end

    return terms
end

function exponent_support(function_::APNFunction)
    n = infer_n(function_)
    terms = expanded_terms(function_)
    if n === nothing
        return Set(term.exponent for term in terms)
    end

    modulus = 2^n - 1
    return Set(normalize_exponent(term.exponent, modulus) for term in terms)
end

function exact_terms(function_::APNFunction, expected::Vector{APNTerm})
    n = infer_n(function_)
    actual_keys = sort([term_key(term, n) for term in expanded_terms(function_)])
    expected_keys = sort([term_key(term, n) for term in expected])
    return actual_keys == expected_keys
end

function trace_has_terms(trace::APNTraceTerm, expected_exponents::Vector{Int})
    actual = sort([term.exponent for term in trace.terms])
    return actual == sort(expected_exponents)
end

function trace_step_matches(trace::APNTraceTerm, step::Int)
    return trace.step === nothing || trace.step == step
end

function has_base_term(function_::APNFunction, exponent::Int)
    n = infer_n(function_)
    return any(term_key(term, n) == term_key(x(exponent), n) for term in function_.terms)
end

function match_c4(function_::APNFunction)
    n = infer_n(function_)
    n === nothing && return nothing

    compact = length(function_.terms) == 1 &&
              has_base_term(function_, 3) &&
              length(function_.traces) == 1 &&
              function_.traces[1].n == n &&
              trace_step_matches(function_.traces[1], 1) &&
              trace_has_terms(function_.traces[1], [9])

    exact_expanded = exact_terms(function_, [x(3); [x(normalize_exponent(9 * 2^j, 2^n - 1)) for j in 0:(n - 1)]])

    (compact || exact_expanded) || return nothing
    return APNFamilyMatch(:C4, Dict(:n => n, :a => 1), true, String[])
end

function match_c5(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && n % 3 == 0) || return nothing

    compact = length(function_.terms) == 1 &&
              has_base_term(function_, 3) &&
              length(function_.traces) == 1 &&
              function_.traces[1].n == n &&
              trace_step_matches(function_.traces[1], 3) &&
              trace_has_terms(function_.traces[1], [9, 18])

    modulus = 2^n - 1
    expected = [x(3)]
    for j in 0:(div(n, 3) - 1)
        push!(expected, x(normalize_exponent(9 * 2^(3j), modulus)))
        push!(expected, x(normalize_exponent(18 * 2^(3j), modulus)))
    end

    exact_expanded = exact_terms(function_, expected)
    (compact || exact_expanded) || return nothing

    return APNFamilyMatch(:C5, Dict(:n => n, :a => 1), true, String[])
end

function match_c6(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && n % 3 == 0) || return nothing

    compact = length(function_.terms) == 1 &&
              has_base_term(function_, 3) &&
              length(function_.traces) == 1 &&
              function_.traces[1].n == n &&
              trace_step_matches(function_.traces[1], 3) &&
              trace_has_terms(function_.traces[1], [18, 36])

    compact || return nothing
    return APNFamilyMatch(:C6, Dict(:n => n, :a => 1), true, String[])
end

function c3_base_exponents(m::Int, i::Int)
    q = 2^m
    n = 2m
    modulus = 2^n - 1
    return Set([
        normalize_exponent(q + 1, modulus),
        normalize_exponent(2^i + 1, modulus),
        normalize_exponent(q * (2^i + 1), modulus),
        normalize_exponent(2^i * q + 1, modulus),
        normalize_exponent(2^i + q, modulus),
    ])
end

function frobenius_orbit(exponent::Int, n::Int)
    modulus = 2^n - 1
    return Set(normalize_exponent(exponent * 2^j, modulus) for j in 0:(n - 1))
end

function support_closure(exponents, n::Int)
    closure = Set{Int}()
    for exponent in exponents
        union!(closure, frobenius_orbit(exponent, n))
    end
    return closure
end

function support_match(family::Symbol, function_::APNFunction, base, parameters;
                       require_base::Bool = true,
                       note::String = "matched by exponent support; coefficient and field-element conditions still need checking")
    n = infer_n(function_)
    n === nothing && return nothing

    support = exponent_support(function_)
    base_support = Set(base)
    support == base_support && return APNFamilyMatch(family, parameters, true, String[])

    closure = support_closure(base_support, n)
    if issubset(support, closure) && (!require_base || issubset(base_support, support))
        return APNFamilyMatch(family, parameters, false, [note])
    end

    return nothing
end

function match_c1_c2(function_::APNFunction)
    n = infer_n(function_)
    n === nothing && return nothing

    for p in (3, 4)
        n % p == 0 || continue
        k = div(n, p)
        gcd(k, 3) == 1 || continue

        for s in 1:(3k - 1)
            gcd(s, 3k) == 1 || continue
            i = mod(s * k, p)
            i == 0 && continue
            m = p - i
            base = [
                normalize_exponent(2^s + 1, 2^n - 1),
                normalize_exponent(2^(i * k) + 2^(m * k + s), 2^n - 1),
            ]
            match = support_match(:C1_C2, function_, base, Dict(:n => n, :p => p, :k => k, :s => s, :i => i, :m => m))
            match !== nothing && return match
        end
    end

    return nothing
end

function match_c3(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && iseven(n)) || return nothing

    m = div(n, 2)
    support = exponent_support(function_)

    for i in 1:(n - 1)
        gcd(i, m) == 1 || continue
        base = c3_base_exponents(m, i)
        support == base && return APNFamilyMatch(:C3, Dict(:n => n, :m => m, :q => 2^m, :i => i), true, String[])

        closure = Set{Int}()
        for exponent in base
            union!(closure, frobenius_orbit(exponent, n))
        end

        if issubset(support, closure) && issubset(base, support)
            return APNFamilyMatch(
                :C3,
                Dict(:n => n, :m => m, :q => 2^m, :i => i),
                false,
                ["support is in the Frobenius-output closure of C3; coefficients/field conditions still need checking"],
            )
        end
    end

    return nothing
end

function match_c7_c8_c9(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && n % 3 == 0) || return nothing

    k = div(n, 3)
    gcd(k, 3) == 1 || return nothing
    modulus = 2^n - 1

    for s in 1:(3k - 1)
        gcd(s, 3k) == 1 || continue
        (k + s) % 3 == 0 || continue

        base = [
            normalize_exponent(2^s + 1, modulus),
            normalize_exponent(2^(n - k) + 2^(k + s), modulus),
            normalize_exponent(2^(n - k) + 1, modulus),
            normalize_exponent(2^s + 2^(k + s), modulus),
        ]
        match = support_match(:C7_C8_C9, function_, base, Dict(:n => n, :k => k, :s => s))
        match !== nothing && return match
    end

    return nothing
end

function match_c10(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && iseven(n)) || return nothing

    m = div(n, 2)
    (m >= 2 && iseven(m)) || return nothing
    modulus = 2^n - 1

    for k in 1:(m - 1)
        gcd(k, m) == 1 || continue
        first_block = [
            2^k + 1,
            2^k + 2^m,
            1 + 2^(m + k),
            2^m * (2^k + 1),
        ]

        for i in 2:2:(n - 1)
            second_block = [exponent * 2^i for exponent in first_block]
            third_block = [2, 1 + 2^m, 2^(m + 1)]
            base = normalize_exponent.([first_block; second_block; third_block], modulus)
            match = support_match(:C10, function_, base, Dict(:n => n, :m => m, :k => k, :i => i), require_base = false)
            match !== nothing && return match
        end
    end

    return nothing
end

function match_c11(function_::APNFunction)
    n = infer_n(function_)
    (n !== nothing && n % 3 == 0) || return nothing

    m = div(n, 3)
    isodd(m) || return nothing
    modulus = 2^n - 1
    base = normalize_exponent.([
        2^(2m + 1) + 1,
        2^(m + 1) + 1,
        2^(2m) + 2,
        2^m + 2,
        3,
    ], modulus)

    return support_match(:C11, function_, base, Dict(:n => n, :m => m), require_base = false)
end

function classify_family(function_::APNFunction)
    matches = APNFamilyMatch[]
    for matcher in (match_c1_c2, match_c3, match_c4, match_c5, match_c6, match_c7_c8_c9, match_c10, match_c11)
        match = matcher(function_)
        match !== nothing && push!(matches, match)
    end
    return matches
end

function classify_family(source::AbstractString; n::Union{Nothing, Int} = nothing)
    return classify_family(parse_apn_function(source; n = n))
end

belongs_to_family(function_::APNFunction, family::Symbol) =
    any(match.family == family for match in classify_family(function_))

belongs_to_family(source::AbstractString, family::Symbol; n::Union{Nothing, Int} = nothing) =
    belongs_to_family(parse_apn_function(source; n = n), family)

function parse_apn_function(source::AbstractString; n::Union{Nothing, Int} = nothing)
    text = replace(source, " " => "")
    isempty(text) && error("Empty APN function")

    parts = split_top_level_sum(text)
    function_n = n
    parsed_parts = Any[]

    for part in parts
        if startswith(part, "Tr")
            trace = parse_trace_term(part)
            function_n === nothing && (function_n = trace.n)
            push!(parsed_parts, trace)
        else
            push!(parsed_parts, parse_apn_term(part))
        end
    end

    return APNFunction(function_n, parsed_parts...)
end

function split_top_level_sum(text::AbstractString)
    parts = String[]
    depth = 0
    start = firstindex(text)

    for index in eachindex(text)
        char = text[index]
        if char == '('
            depth += 1
        elseif char == ')'
            depth -= 1
            depth < 0 && error("Unbalanced parentheses in $text")
        elseif char == '+' && depth == 0
            push!(parts, String(text[start:prevind(text, index)]))
            start = nextind(text, index)
        end
    end

    depth == 0 || error("Unbalanced parentheses in $text")
    push!(parts, String(text[start:lastindex(text)]))
    return parts
end

function parse_trace_term(text::AbstractString)
    open_index = findfirst(==('('), text)
    close_index = findlast(==(')'), text)
    (open_index !== nothing && close_index == lastindex(text)) || error("Invalid trace term: $text")

    header = text[firstindex(text):prevind(text, open_index)]
    body = text[nextind(text, open_index):prevind(text, close_index)]
    header_match = match(r"^Tr([0-9]+)(?:\^([0-9]+))?$", header)
    header_match !== nothing || error("Invalid trace header: $header")

    n = parse(Int, header_match.captures[1])
    step = header_match.captures[2] === nothing ? nothing : parse(Int, header_match.captures[2])
    terms = [parse_apn_term(part) for part in split_top_level_sum(body)]
    return APNTraceTerm(n, step, terms)
end

function parse_apn_term(text::AbstractString)
    term_match = match(r"^(?:(p)([0-9]*))?x([0-9]+)$", text)
    term_match !== nothing || error("Invalid APN term: $text")

    if term_match.captures[1] === nothing
        coefficient = ONE_COEFF
    else
        exponent_text = term_match.captures[2]
        coefficient_power = isempty(exponent_text) ? 1 : parse(Int, exponent_text)
        coefficient = PowerCoefficient(:p, coefficient_power)
    end

    exponent = parse(Int, term_match.captures[3])
    return APNTerm(coefficient, exponent)
end

function family_c1_c2()
    error("C1-C2 recognition is not implemented yet")
end

function family_c3(n::Int, i::Int)
    iseven(n) || error("C3 requires n = 2m")
    m = div(n, 2)
    gcd(i, m) == 1 || error("C3 requires gcd(i, m) = 1")
    return c3_base_exponents(m, i)
end

function family_c4(n::Int)
    return APNFunction(n, x(3), Tr(n, x(9); step = 1))
end

function family_c5(n::Int)
    n % 3 == 0 || error("C5 requires 3 | n")
    return APNFunction(n, x(3), Tr(n, x(9), x(18); step = 3))
end

function family_c6(n::Int)
    n % 3 == 0 || error("C6 requires 3 | n")
    return APNFunction(n, x(3), Tr(n, x(18), x(36); step = 3))
end

function family_c7_c8_c9()
    error("C7-C9 recognition is not implemented yet")
end

function family_c10()
    error("C10 recognition is not implemented yet")
end

function family_c11()
    error("C11 recognition is not implemented yet")
end

function family_gold(n::Int, i::Int, d::Int)
    gcd(i, n) == 1 || error("i and n must be coprime")
    d == 2^i + 1 || error("d must be equal to 2^i + 1")
    return true
end

function family_kasami(n::Int, i::Int, d::Int)
    gcd(i, n) == 1 || error("i and n must be coprime")
    d == 2^(2i) - 2^i + 1 || error("d must be equal to 2^(2i) - 2^i + 1")
    return true
end

function family_welsh()
    error("Welch recognition is not implemented yet")
end

function family_niho()
    error("Niho recognition is not implemented yet")
end

function family_inverse()
    error("Inverse recognition is not implemented yet")
end

function family_dobbertin()
    error("Dobbertin recognition is not implemented yet")
end
