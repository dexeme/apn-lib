const PRECOMPUTED_TUPLE_CACHE = Dict{Int, Vector{Vector{Int}}}()
const PRECOMPUTED_TUPLE_MATRIX_CACHE = Dict{Int, Any}()

default_tuples_dir() = joinpath(@__DIR__, "..", "..", "tuples")

function parse_c_tuple_rows(filename::String)
    text = read(filename, String)
    body_match = match(r"int\s+AllTuples\s*\[[^\]]+\]\s*\[[^\]]+\]\s*=\s*\{(.*)\};"s, text)
    body_match === nothing && error("Could not find AllTuples array in $filename")

    rows = Vector{Vector{Int}}()

    for row_match in eachmatch(r"\{([^{}]*)\}", body_match.captures[1])
        values = Int[]
        for value_match in eachmatch(r"-?\d+", row_match.captures[1])
            push!(values, parse(Int, value_match.match))
        end
        push!(rows, values)
    end

    return rows
end

function tuple_constants_name(n::Int, index::Int)
    return "ALL_TUPLES_$(n)_$(index)"
end

function write_julia_tuple_constants(file, n::Int, rows::Vector{Vector{Int}})
    write(file, "# This file is generated from tuples/AllTuples$n.h.\n")
    write(file, "# Regenerate it with generate_tuple_constants_file($n).\n\n")
    write(file, "const N_TUPLES_$n = $(length(rows))\n\n")

    names = String[]
    for index in 1:length(rows)
        name = tuple_constants_name(n, index)
        push!(names, name)
        write(file, "const $name = Int[$(join(rows[index], ", "))]\n\n")
    end

    write(file, "const ALL_TUPLES_$n = Vector{Int}[$(join(names, ", "))]\n")
end

function generate_tuple_constants_file(n::Int; tuples_dir::String = default_tuples_dir())
    source_filename = joinpath(tuples_dir, "AllTuples$n.h")
    output_filename = joinpath(tuples_dir, "AllTuples$n.jl")

    rows = parse_c_tuple_rows(source_filename)

    open(output_filename, "w") do file
        write_julia_tuple_constants(file, n, rows)
    end

    return output_filename
end

function generate_tuple_constants_files(ns; tuples_dir::String = default_tuples_dir())
    return [generate_tuple_constants_file(n, tuples_dir = tuples_dir) for n in ns]
end

function include_constants_namespace(namespace_name::Symbol, filename::String)
    namespace = Module(namespace_name)
    Base.include(namespace, filename)
    return namespace
end

function get_constant(namespace::Module, name::Symbol)
    return Base.invokelatest(getfield, namespace, name)
end

function load_precomputed_tuple_constants(n::Int; tuples_dir::String = default_tuples_dir())
    if haskey(PRECOMPUTED_TUPLE_CACHE, n)
        return PRECOMPUTED_TUPLE_CACHE[n]
    end

    filename = joinpath(tuples_dir, "AllTuples$n.jl")
    if !isfile(filename)
        matrix_filename = joinpath(tuples_dir, "AllTuplesMatrices$n.jl")
        if isfile(matrix_filename)
            rows = [block.tuple for block in load_precomputed_tuple_matrix_constants(n, tuples_dir = tuples_dir)]
            PRECOMPUTED_TUPLE_CACHE[n] = rows
            return rows
        end

        error("Precomputed constants file not found: $filename")
    end

    namespace = include_constants_namespace(Symbol("APNLibPrecomputedTuples", n), filename)

    rows = get_constant(namespace, Symbol("ALL_TUPLES_$n"))
    PRECOMPUTED_TUPLE_CACHE[n] = rows

    return rows
end

function precomputed_tuple_row(n::Int, index::Int; tuples_dir::String = default_tuples_dir())
    rows = load_precomputed_tuple_constants(n, tuples_dir = tuples_dir)
    1 <= index <= length(rows) || error("Tuple index must be between 1 and $(length(rows)) for n = $n")

    return rows[index]
end

function precomputed_tuple_sboxes(n::Int, index::Int; tuples_dir::String = default_tuples_dir())
    tuple_lut = precomputed_tuple_row(n, index, tuples_dir = tuples_dir)
    sbox_size = 2^n
    length(tuple_lut) == 2 * sbox_size || error("Tuple row for n = $n must have $(2 * sbox_size) entries")

    lut_B = tuple_lut[1:sbox_size]
    lut_A = tuple_lut[(sbox_size + 1):end]

    return lut_A, lut_B
end

function lut_to_matrix(lut::AbstractVector{<:Integer}, n::Int)::Matrix{Int}
    length(lut) == 2^n || error("LUT for n = $n must have $(2^n) entries")

    matrix = zeros(Int, n, n)

    for col in 1:n
        basis_vector = 2^(col - 1)
        value = lut[basis_vector + 1]
        matrix[:, col] = digits(value, base = 2, pad = n)[1:n]
    end

    return matrix
end

function extract_matrices_from_tuple_lut(tuple_lut::AbstractVector{<:Integer}, n::Int)
    sbox_size = 2^n
    length(tuple_lut) == 2 * sbox_size || error("Tuple LUT for n = $n must have $(2 * sbox_size) entries")

    lut_B = tuple_lut[1:sbox_size]
    lut_A = tuple_lut[(sbox_size + 1):end]

    A = lut_to_matrix(lut_A, n)
    B = lut_to_matrix(lut_B, n)

    return A, B
end

function extract_matrices(tuple_lut::AbstractVector{<:Integer}, n::Int)
    return extract_matrices_from_tuple_lut(tuple_lut, n)
end

function matrix_literal(matrix::AbstractMatrix{<:Integer})
    rows = String[]

    for row in 1:size(matrix, 1)
        push!(rows, join((string(Int(matrix[row, col])) for col in 1:size(matrix, 2)), " "))
    end

    return "Int[" * join(rows, "; ") * "]"
end

function int_vector_literal(values::AbstractVector{<:Integer})
    return "Int[$(join((string(Int(value)) for value in values), ", "))]"
end

function write_julia_tuple_matrix_constants(file, n::Int, rows::Vector{Vector{Int}})
    write(file, "# This file is generated from tuples/AllTuples$n.h.\n")
    write(file, "# Regenerate it with generate_tuple_matrix_constants_file($n).\n")
    write(file, "# Each class is stored as: tuple LUT, matrix A, matrix B.\n\n")
    write(file, "const N_TUPLE_MATRIX_BLOCKS_$n = $(length(rows))\n\n")

    for index in 1:length(rows)
        tuple_name = "ALL_TUPLES_$(n)_$(index)_TUPLE"
        a_name = "ALL_TUPLES_$(n)_$(index)_A"
        b_name = "ALL_TUPLES_$(n)_$(index)_B"
        A, B = extract_matrices_from_tuple_lut(rows[index], n)

        write(file, "const $tuple_name = Int[$(join(rows[index], ", "))]\n")
        write(file, "const $a_name = $(matrix_literal(A))\n")
        write(file, "const $b_name = $(matrix_literal(B))\n")
        write(file, "\n")
    end
end

function generate_tuple_matrix_constants_file(n::Int; tuples_dir::String = default_tuples_dir())
    source_filename = joinpath(tuples_dir, "AllTuples$n.h")
    output_filename = joinpath(tuples_dir, "AllTuplesMatrices$n.jl")

    rows = parse_c_tuple_rows(source_filename)

    open(output_filename, "w") do file
        write_julia_tuple_matrix_constants(file, n, rows)
    end

    return output_filename
end

function generate_tuple_matrix_constants_files(ns; tuples_dir::String = default_tuples_dir())
    return [generate_tuple_matrix_constants_file(n, tuples_dir = tuples_dir) for n in ns]
end

function tuple_matrix_block_from_namespace(namespace::Module, n::Int, index::Int)
    search_name = Symbol("ALL_TUPLES_$(n)_$(index)_SEARCH")
    search_result = isdefined(namespace, search_name) ? get_constant(namespace, search_name) : nothing

    return (
        tuple = get_constant(namespace, Symbol("ALL_TUPLES_$(n)_$(index)_TUPLE")),
        A = get_constant(namespace, Symbol("ALL_TUPLES_$(n)_$(index)_A")),
        B = get_constant(namespace, Symbol("ALL_TUPLES_$(n)_$(index)_B")),
        search = search_result,
    )
end

function load_precomputed_tuple_matrix_constants(n::Int; tuples_dir::String = default_tuples_dir())
    if haskey(PRECOMPUTED_TUPLE_MATRIX_CACHE, n)
        return PRECOMPUTED_TUPLE_MATRIX_CACHE[n]
    end

    filename = joinpath(tuples_dir, "AllTuplesMatrices$n.jl")
    isfile(filename) || error("Precomputed matrix constants file not found: $filename")

    namespace = include_constants_namespace(Symbol("APNLibPrecomputedTupleMatrices", n), filename)

    n_blocks = get_constant(namespace, Symbol("N_TUPLE_MATRIX_BLOCKS_$n"))
    blocks = [tuple_matrix_block_from_namespace(namespace, n, index) for index in 1:n_blocks]
    PRECOMPUTED_TUPLE_MATRIX_CACHE[n] = blocks

    return blocks
end

function precomputed_tuple_matrices(n::Int, index::Int; tuples_dir::String = default_tuples_dir())
    blocks = load_precomputed_tuple_matrix_constants(n, tuples_dir = tuples_dir)
    1 <= index <= length(blocks) || error("Tuple index must be between 1 and $(length(blocks)) for n = $n")

    block = blocks[index]
    return block.A, block.B
end

function search_result_constant_name(n::Int, class_index::Int)::String
    return "ALL_TUPLES_$(n)_$(class_index)_SEARCH"
end

function tuple_matrix_constants_filename(n::Int)
    return joinpath(default_tuples_dir(), "AllTuplesMatrices$n.jl")
end

function save_search_result_constant(sbox::AbstractVector{<:Integer}, n::Int, class_index::Int)
    space_size = 2^n
    length(sbox) == space_size || error("Search result for n = $n must have $space_size entries")
    1 <= class_index || error("class_index must be positive")

    filename = tuple_matrix_constants_filename(n)
    isfile(filename) || error("Precomputed matrix constants file not found: $filename")

    constant_name = search_result_constant_name(n, class_index)
    constant_line = "const $constant_name = $(int_vector_literal(sbox))"
    text = read(filename, String)

    existing_pattern = Regex("(?m)^const\\s+$constant_name\\s*=\\s*Int\\[[^\\n]*\\]\$")
    if occursin(existing_pattern, text)
        text = replace(text, existing_pattern => constant_line)
    else
        b_name = "ALL_TUPLES_$(n)_$(class_index)_B"
        b_pattern = Regex("(?m)^const\\s+$b_name\\s*=\\s*Int\\[[^\\n]*\\]\\n")

        if occursin(b_pattern, text)
            text = replace(text, b_pattern => s -> s * constant_line * "\n")
        else
            text *= "\n$constant_line\n"
        end
    end

    write(filename, text)
    delete!(PRECOMPUTED_TUPLE_MATRIX_CACHE, n)

    return filename
end

function tuple_to_sbox_row(tuple)
    A, B = tuple
    return vcat(matrix_to_sbox(A), matrix_to_sbox(B))
end

function tuples_to_sbox_rows(tuples)
    return [tuple_to_sbox_row(tuple) for tuple in tuples]
end

function write_c_array(file, rows)
    if isempty(rows)
        write(file, "int AllTuples[0][0] = {};\n")
        return
    end

    row_size = length(rows[1])

    write(file, "#define N_TUPLES $(length(rows))\n\n")
    write(file, "int AllTuples[N_TUPLES][$row_size] = {\n")

    for i in 1:length(rows)
        row = rows[i]
        row_text = join(row, ", ")

        if i < length(rows)
            write(file, "    {$row_text},\n")
        else
            write(file, "    {$row_text}\n")
        end
    end

    write(file, "};\n")
end

function generate_tuples_file(tuples, filename::String = "tuples/AllTuples.h")
    rows = tuples_to_sbox_rows(tuples)

    open(filename, "w") do file
        write_c_array(file, rows)
    end
end
