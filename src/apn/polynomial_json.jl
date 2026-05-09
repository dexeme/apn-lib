using JSON3
using Nemo

function parse_binary_field_degree(field::AbstractString)
    field_match = match(r"^GF\(2\^([0-9]+)\)$", replace(field, " " => ""))
    field_match !== nothing || error("Only fields written as GF(2^n) are supported")
    return parse(Int, field_match.captures[1])
end

function parse_binary_modulus(modulus::AbstractString)
    compact = replace(modulus, " " => "")
    variable_match = match(r"^([A-Za-z_][A-Za-z_0-9]*)", compact)
    variable_match !== nothing || error("Could not infer modulus variable from $modulus")
    variable_name = variable_match.captures[1]

    base_field = GF(2)
    polynomial_ring_, variable = polynomial_ring(base_field, variable_name)
    polynomial = zero(polynomial_ring_)

    for raw_term in split(compact, "+")
        if raw_term == "1"
            polynomial += one(polynomial_ring_)
        elseif raw_term == variable_name
            polynomial += variable
        else
            term_match = match(Regex("^" * variable_name * "\\^([0-9]+)\$"), raw_term)
            term_match !== nothing || error("Invalid modulus term: $raw_term")
            polynomial += variable^parse(Int, term_match.captures[1])
        end
    end

    return polynomial, variable_name
end

function polynomial_json_data(json_text::AbstractString)
    return JSON3.read(json_text)
end

function json_property(data, name::Symbol)
    return getproperty(data, name)
end

function json_property(data::AbstractDict, name::Symbol)
    return data[String(name)]
end

function json_has_property(data, name::Symbol)
    return hasproperty(data, name)
end

function json_has_property(data::AbstractDict, name::Symbol)
    return haskey(data, String(name))
end

function polynomial_json_terms(data)
    terms = json_property(data, :terms)
    return collect(terms)
end

function polynomial_json_dimension(data)
    return parse_binary_field_degree(String(json_property(data, :field)))
end

function canonical_polynomial_json(data)
    return JSON3.write(data)
end

function parse_polynomial_expression_terms(expression::AbstractString)
    terms = Dict{String, Int}[]
    compact_terms = split(replace(expression, " " => ""), "+")

    for raw_term in compact_terms
        isempty(raw_term) && continue
        term_match = match(r"^(?:g([0-9]*))?x([0-9]+)$", raw_term)
        term_match !== nothing || error("Invalid polynomial term: $raw_term")

        coefficient_text = term_match.captures[1]
        exponent = parse(Int, term_match.captures[2])

        if coefficient_text === nothing
            push!(terms, Dict("exp" => exponent, "coef" => 1))
        else
            coefficient_power = isempty(coefficient_text) ? 1 : parse(Int, coefficient_text)
            push!(terms, Dict("exp" => exponent, "coef_power" => coefficient_power))
        end
    end

    isempty(terms) && error("Polynomial expression has no terms: $expression")
    return terms
end

function polynomial_expression_json(expression::AbstractString;
                                    dimension::Integer,
                                    modulus::AbstractString,
                                    basis::AbstractString = "power")
    return Dict(
        "field" => "GF(2^$(Int(dimension)))",
        "basis" => String(basis),
        "modulus" => String(modulus),
        "terms" => parse_polynomial_expression_terms(expression),
    )
end

function build_polynomial_from_json(json_text::AbstractString)
    return build_polynomial_from_json(polynomial_json_data(json_text))
end

function build_polynomial_from_json(data)
    field_text = String(json_property(data, :field))
    basis = String(json_property(data, :basis))
    modulus_text = String(json_property(data, :modulus))
    terms = polynomial_json_terms(data)

    basis == "power" || error("Only basis = \"power\" is supported")

    n = parse_binary_field_degree(field_text)
    modulus, field_variable_name = parse_binary_modulus(modulus_text)
    degree(modulus) == n || error("Modulus degree $(degree(modulus)) does not match field degree $n")

    field = GF(modulus, field_variable_name)
    polynomial_ring_, x_variable = polynomial_ring(field, "x")
    polynomial = zero(polynomial_ring_)

    for term in terms
        exponent = Int(json_property(term, :exp))
        coefficient = if json_has_property(term, :coef_power)
            gen(field)^Int(json_property(term, :coef_power))
        else
            int_to_field_element(Int(json_property(term, :coef)), field, n)
        end
        polynomial += coefficient * x_variable^exponent
    end

    return polynomial
end

function symbolic_apn_function_from_json(data)
    n = polynomial_json_dimension(data)
    terms = [x(Int(json_property(term, :exp))) for term in polynomial_json_terms(data)]
    return APNFunction(n, terms...)
end

function normalized_support_from_json(data)
    n = polynomial_json_dimension(data)
    modulus = 2^n - 1
    support = sort(unique(normalize_exponent(Int(json_property(term, :exp)), modulus) for term in polynomial_json_terms(data)))
    return JSON3.write(support)
end
