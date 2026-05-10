function sqlite_identifier(name::AbstractString)
    occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", name) || error("Invalid SQL identifier: $name")
    return name
end

function upsert_table_values!(db::SQLite.DB,
                              table::AbstractString,
                              key_column::AbstractString,
                              key_value,
                              values)
    value_pairs = collect(values)
    isempty(value_pairs) && return nothing

    table_name = sqlite_identifier(table)
    key_name = sqlite_identifier(key_column)
    columns = [sqlite_identifier(String(column)) for (column, _) in value_pairs]
    placeholders = join(fill("?", length(columns) + 1), ", ")
    insert_columns = join([key_name; columns], ", ")
    updates = join(["$column = excluded.$column" for column in columns], ", ")
    value_tuple = Tuple(value for (_, value) in value_pairs)

    DBInterface.execute(
        db,
        """
        INSERT INTO $table_name ($insert_columns)
        VALUES ($placeholders)
        ON CONFLICT($key_name) DO UPDATE SET $updates
        """,
        (key_value, value_tuple...),
    )

    return key_value
end

function update_table_values!(db::SQLite.DB,
                              table::AbstractString,
                              values;
                              where_column::AbstractString,
                              where_value)
    value_pairs = collect(values)
    isempty(value_pairs) && return nothing

    table_name = sqlite_identifier(table)
    where_name = sqlite_identifier(where_column)
    assignments = join(["$(sqlite_identifier(String(column))) = ?" for (column, _) in value_pairs], ", ")
    value_tuple = Tuple(value for (_, value) in value_pairs)

    DBInterface.execute(
        db,
        "UPDATE $table_name SET $assignments WHERE $where_name = ?",
        (value_tuple..., where_value),
    )

    return where_value
end

function json_table(table_json::AbstractString)
    rows = JSON3.read(table_json)
    length(rows) >= 1 || error("JSON table must have a header row")
    return String.(collect(rows[1])), [collect(row) for row in rows[2:end]]
end

function json_row_dict(header, row)
    length(row) == length(header) || error("JSON table row has $(length(row)) values, expected $(length(header))")
    return Dict(String(column) => value for (column, value) in zip(header, row))
end

function split_int_values(value)
    return [parse(Int, strip(part)) for part in split(string(value), ",") if !isempty(strip(part))]
end

function require_json_columns(header, columns)
    missing = setdiff(String.(collect(columns)), String.(header))
    isempty(missing) || error("Missing JSON table columns: $(join(missing, ", "))")
    return nothing
end
