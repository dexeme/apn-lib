using Nemo

function _coefficient_to_field_element(::OneCoefficient, field::FqField, n::Int)
    return one(field)
end

function _coefficient_to_field_element(coefficient::PowerCoefficient, field::FqField, n::Int)
    check_binary_extension_field(field, n)

    modulus = space_size(n) - 1
    modulus > 0 || return one(field)

    generator = gen(field)
    exponent = mod(coefficient.exponent, modulus)

    return generator^exponent
end

function _evaluate_term(term::APNTerm, x_value::FqFieldElem, field::FqField, n::Int)
    coefficient = _coefficient_to_field_element(term.coefficient, field, n)
    return coefficient * x_value^term.exponent
end

function _relative_trace_to_field(element::FqFieldElem, step::Int)::FqFieldElem
    field = parent(element)
    n = degree(field)

    step > 0 || throw(ArgumentError("relative trace step must be positive"))
    n % step == 0 || throw(ArgumentError("relative trace step must divide the field degree"))

    result = zero(field)
    current = element

    for _ in 1:div(n, step)
        result += current
        current = current^(1 << step)
    end

    return result
end

function _evaluate_trace_term(trace::APNTraceTerm, x_value::FqFieldElem, field::FqField, n::Int)
    trace.n == n ||
        throw(ArgumentError("trace dimension $(trace.n) does not match function dimension $n"))

    value = zero(field)

    for term in trace.terms
        value += _evaluate_term(term, x_value, field, n)
    end

    if trace.step === nothing
        return absolute_trace_to_field(value)
    end

    return _relative_trace_to_field(value, trace.step)
end

function evaluate(function_::APNFunction, x_value::FqFieldElem, field::FqField)
    function_.n !== nothing ||
        throw(ArgumentError("cannot evaluate APNFunction without dimension n"))
    isempty(function_.components) ||
        throw(ArgumentError("cannot evaluate APNFunction with catalogue components without an evaluation context"))

    n = function_.n
    check_binary_extension_field(field, n)

    parent(x_value) == field ||
        throw(ArgumentError("x_value must belong to the provided field"))

    result = zero(field)

    for term in function_.terms
        result += _evaluate_term(term, x_value, field, n)
    end

    for trace in function_.traces
        result += _evaluate_trace_term(trace, x_value, field, n)
    end

    return result
end

function apn_to_lut(function_::APNFunction)::Vector{Int}
    function_.n !== nothing ||
        throw(ArgumentError("cannot convert APNFunction without dimension n"))

    n = function_.n
    field = GF(2, n, "g")
    lookup = field_element_lookup(field, n)

    lut = Vector{Int}(undef, space_size(n))

    for x_int in 0:(space_size(n) - 1)
        x_value = int_to_field_element(x_int, field, n)
        y_value = evaluate(function_, x_value, field)

        haskey(lookup, y_value) ||
            throw(ArgumentError("function output does not belong to GF(2^$n)"))

        lut[x_int + 1] = lookup[y_value]
    end

    return lut
end

function apn_with_dimension(function_::APNFunction, n::Int)::APNFunction
    if function_.n === nothing
        return APNFunction(n, function_.id, function_)
    end

    function_.n == n ||
        throw(ArgumentError("function dimension $(function_.n) does not match n = $n"))

    return function_
end

apn_to_lut(function_::APNFunction, n::Int)::Vector{Int} =
    apn_to_lut(apn_with_dimension(function_, n))

univariate_to_lut(function_::APNFunction, n::Int)::Vector{Int} =
    apn_to_lut(function_, n)

univariate_to_graph(function_::APNFunction, n::Int)::Vector{Tuple{Int, Int}} =
    lut_to_graph(univariate_to_lut(function_, n), n)

univariate_to_anf(function_::APNFunction, n::Int)::ANFVector =
    lut_to_anf(univariate_to_lut(function_, n), n)
