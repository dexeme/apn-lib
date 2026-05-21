function pretty_row_value(row, column::Symbol)
    if column == :delta && (:found in propertynames(row)) && (:expected in propertynames(row))
        return getproperty(row, :found) - getproperty(row, :expected)
    end

    return getproperty(row, column)
end

function pretty_cell(value; digits::Int = 6)
    if value isa AbstractFloat
        return string(round(value; digits = digits))
    end

    if value isa Bool
        return value ? "yes" : "no"
    end

    return string(value)
end

function pretty_header(column::Symbol, headers)
    return get(headers, column, uppercasefirst(String(column)))
end

function pretty_alignment(column::Symbol, alignments)
    alignment = get(alignments, column, :left)

    if alignment == :right
        return "---:"
    elseif alignment == :center
        return ":---:"
    elseif alignment == :left
        return ":---"
    end

    error("unsupported table alignment: $alignment")
end

function pretty_print_table(rows;
                            columns::Vector{Symbol},
                            headers = Dict{Symbol, String}(),
                            alignments = Dict{Symbol, Symbol}(),
                            io::IO = stdout,
                            digits::Int = 6)
    println(io, "| ", join((pretty_header(column, headers) for column in columns), " | "), " |")
    println(io, "| ", join((pretty_alignment(column, alignments) for column in columns), " | "), " |")

    for row in rows
        values = [pretty_cell(pretty_row_value(row, column); digits = digits) for column in columns]
        println(io, "| ", join(values, " | "), " |")
        flush(io)
    end

    return nothing
end

function pretty_print_mismatches(rows;
                                 columns::Vector{Symbol} = [:n, :id, :equation, :found, :expected, :delta],
                                 headers = Dict{Symbol, String}(:id => "ID"),
                                 alignments = Dict{Symbol, Symbol}(:n => :right,
                                                                   :found => :right,
                                                                   :expected => :right,
                                                                   :delta => :right),
                                 ok_column::Symbol = :ok,
                                 io::IO = stdout,
                                 digits::Int = 6)
    mismatches = [row for row in rows if (ok_column in propertynames(row)) && !getproperty(row, ok_column)]
    isempty(mismatches) && return nothing

    println(io)
    println(io, "Mismatches:")
    pretty_print_table(mismatches;
                       columns = columns,
                       headers = headers,
                       alignments = alignments,
                       io = io,
                       digits = digits)
    return nothing
end

function pretty_print_experiment_result(result::ExperimentResult;
                                        columns::Vector{Symbol} = result.spec.columns,
                                        headers = Dict{Symbol, String}(),
                                        alignments = Dict{Symbol, Symbol}(),
                                        mismatch_columns::Union{Nothing, Vector{Symbol}} = nothing,
                                        io::IO = stdout,
                                        digits::Int = 6)
    pretty_print_table(result.rows;
                       columns = columns,
                       headers = headers,
                       alignments = alignments,
                       io = io,
                       digits = digits)

    if mismatch_columns !== nothing
        pretty_print_mismatches(result.rows;
                                columns = mismatch_columns,
                                headers = headers,
                                alignments = alignments,
                                io = io,
                                digits = digits)
    end

    return nothing
end
