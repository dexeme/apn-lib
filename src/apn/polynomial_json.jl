using Nemo

struct JsonCursor
    text::String
    index::Int
end

function skip_json_space(cursor::JsonCursor)
    index = cursor.index
    while index <= lastindex(cursor.text) && cursor.text[index] in (' ', '\n', '\r', '\t')
        index = nextind(cursor.text, index)
    end
    return JsonCursor(cursor.text, index)
end

function expect_json_char(cursor::JsonCursor, expected::Char)
    cursor = skip_json_space(cursor)
    cursor.index <= lastindex(cursor.text) || error("Expected '$expected', got end of JSON")
    cursor.text[cursor.index] == expected || error("Expected '$expected', got '$(cursor.text[cursor.index])'")
    return JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
end

function parse_json_string(cursor::JsonCursor)
    cursor = expect_json_char(cursor, '"')
    buffer = IOBuffer()
    index = cursor.index

    while index <= lastindex(cursor.text)
        char = cursor.text[index]
        if char == '"'
            return String(take!(buffer)), JsonCursor(cursor.text, nextind(cursor.text, index))
        elseif char == '\\'
            escape_index = nextind(cursor.text, index)
            escape_index <= lastindex(cursor.text) || error("Invalid JSON string escape")
            escaped = cursor.text[escape_index]
            if escaped == '"' || escaped == '\\' || escaped == '/'
                print(buffer, escaped)
            elseif escaped == 'n'
                print(buffer, '\n')
            elseif escaped == 'r'
                print(buffer, '\r')
            elseif escaped == 't'
                print(buffer, '\t')
            else
                error("Unsupported JSON escape: \\$escaped")
            end
            index = nextind(cursor.text, escape_index)
        else
            print(buffer, char)
            index = nextind(cursor.text, index)
        end
    end

    error("Unterminated JSON string")
end

function parse_json_number(cursor::JsonCursor)
    cursor = skip_json_space(cursor)
    start = cursor.index
    index = cursor.index

    if index <= lastindex(cursor.text) && cursor.text[index] == '-'
        index = nextind(cursor.text, index)
    end

    while index <= lastindex(cursor.text) && isdigit(cursor.text[index])
        index = nextind(cursor.text, index)
    end

    start == index && error("Expected JSON integer")
    return parse(Int, cursor.text[start:prevind(cursor.text, index)]), JsonCursor(cursor.text, index)
end

function parse_json_array(cursor::JsonCursor)
    cursor = expect_json_char(cursor, '[')
    values = Any[]
    cursor = skip_json_space(cursor)

    if cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == ']'
        return values, JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
    end

    while true
        value, cursor = parse_json_value(cursor)
        push!(values, value)
        cursor = skip_json_space(cursor)

        if cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == ','
            cursor = JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
        elseif cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == ']'
            return values, JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
        else
            error("Expected ',' or ']' in JSON array")
        end
    end
end

function parse_json_object(cursor::JsonCursor)
    cursor = expect_json_char(cursor, '{')
    object = Dict{String, Any}()
    cursor = skip_json_space(cursor)

    if cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == '}'
        return object, JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
    end

    while true
        key, cursor = parse_json_string(cursor)
        cursor = expect_json_char(cursor, ':')
        value, cursor = parse_json_value(cursor)
        object[key] = value
        cursor = skip_json_space(cursor)

        if cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == ','
            cursor = JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
        elseif cursor.index <= lastindex(cursor.text) && cursor.text[cursor.index] == '}'
            return object, JsonCursor(cursor.text, nextind(cursor.text, cursor.index))
        else
            error("Expected ',' or '}' in JSON object")
        end
    end
end

function parse_json_value(cursor::JsonCursor)
    cursor = skip_json_space(cursor)
    cursor.index <= lastindex(cursor.text) || error("Unexpected end of JSON")
    char = cursor.text[cursor.index]

    char == '{' && return parse_json_object(cursor)
    char == '[' && return parse_json_array(cursor)
    char == '"' && return parse_json_string(cursor)
    (char == '-' || isdigit(char)) && return parse_json_number(cursor)

    error("Unsupported JSON value starting with '$char'")
end

function parse_json_object(text::AbstractString)
    value, cursor = parse_json_value(JsonCursor(String(text), firstindex(String(text))))
    cursor = skip_json_space(cursor)
    cursor.index > lastindex(cursor.text) || error("Unexpected trailing JSON content")
    value isa Dict{String, Any} || error("Top-level JSON value must be an object")
    return value
end

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

function build_polynomial_from_json(json_text::AbstractString)
    data = parse_json_object(json_text)

    field_text = get(data, "field", nothing)
    basis = get(data, "basis", nothing)
    modulus_text = get(data, "modulus", nothing)
    terms = get(data, "terms", nothing)

    field_text isa AbstractString || error("JSON field must be a string, for example GF(2^7)")
    basis == "power" || error("Only basis = \"power\" is supported")
    modulus_text isa AbstractString || error("JSON modulus must be a string")
    terms isa Vector || error("JSON terms must be an array")

    n = parse_binary_field_degree(field_text)
    modulus, field_variable_name = parse_binary_modulus(modulus_text)
    degree(modulus) == n || error("Modulus degree $(degree(modulus)) does not match field degree $n")

    field = GF(modulus, field_variable_name)
    polynomial_ring_, x_variable = polynomial_ring(field, "x")
    polynomial = zero(polynomial_ring_)

    for term in terms
        term isa Dict{String, Any} || error("Each term must be a JSON object")
        exponent = get(term, "exp", nothing)
        coefficient = get(term, "coef", nothing)
        exponent isa Integer || error("Term exp must be an integer")
        coefficient isa Integer || error("Term coef must be an integer")
        polynomial += int_to_field_element(coefficient, field, n) * x_variable^exponent
    end

    return polynomial
end
