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
    return any(term_key(term, n) == term_key(monomial_expr(exponent), n) for term in function_.terms)
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
