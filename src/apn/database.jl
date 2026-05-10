using DBInterface
using JSON3
using SQLite

include(joinpath(@__DIR__, "..", "..", "db", "config.jl"))

const DEFAULT_APN_DB_PATH = APN_DB_PATH
const DEFAULT_APN_SCHEMA_PATH = APN_SCHEMA_PATH

function open_apn_database(path::AbstractString = DEFAULT_APN_DB_PATH)
    return SQLite.DB(path)
end

function init_apn_database!(path::AbstractString = DEFAULT_APN_DB_PATH;
                            schema_path::AbstractString = DEFAULT_APN_SCHEMA_PATH)
    mkpath(dirname(path))
    db = SQLite.DB(path)
    try
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        execute_schema!(db, read(schema_path, String))
        return db
    catch
        SQLite.close(db)
        rethrow()
    end
end

function with_apn_database(f::Function, path::AbstractString = DEFAULT_APN_DB_PATH)
    db = init_apn_database!(path)
    try
        return f(db)
    finally
        SQLite.close(db)
    end
end

function execute_schema!(db::SQLite.DB, schema::AbstractString)
    for statement in split(schema, ";")
        sql = strip(statement)
        isempty(sql) && continue
        DBInterface.execute(db, sql)
    end
    migrate_apn_function_columns!(db)
end

function sqlite_column_exists(db::SQLite.DB, table::AbstractString, column::AbstractString)
    for row in DBInterface.execute(db, "PRAGMA table_info($table)")
        String(row.name) == column && return true
    end
    return false
end

function sqlite_index_exists(db::SQLite.DB, index_name::AbstractString)
    existing = sqlite_first_value(
        db,
        "SELECT 1 FROM sqlite_master WHERE type = 'index' AND name = ?",
        (index_name,),
    )
    return existing !== nothing
end

function migrate_apn_function_columns!(db::SQLite.DB)
    if !sqlite_column_exists(db, "apn_function", "local_id")
        DBInterface.execute(db, "ALTER TABLE apn_function ADD COLUMN local_id INTEGER")
        DBInterface.execute(db, "UPDATE apn_function SET local_id = id WHERE local_id IS NULL")
    end

    if !sqlite_column_exists(db, "apn_function", "equivalent_to")
        DBInterface.execute(db, "ALTER TABLE apn_function ADD COLUMN equivalent_to TEXT")
    end

    if !sqlite_column_exists(db, "apn_function", "walsh_spectrum")
        DBInterface.execute(db, "ALTER TABLE apn_function ADD COLUMN walsh_spectrum TEXT")
    end

    if !sqlite_index_exists(db, "idx_apn_function_dimension_local_id")
        DBInterface.execute(
            db,
            """
            CREATE UNIQUE INDEX IF NOT EXISTS idx_apn_function_dimension_local_id
            ON apn_function (dimension, local_id)
            WHERE local_id IS NOT NULL
            """,
        )
    end
end

function apn_family_db_name(family::Symbol)
    family == :C1_C2 && return "C1-C2"
    family == :C7_C8_C9 && return "C7-C9"
    return String(family)
end

function sqlite_first_value(db::SQLite.DB, sql::AbstractString, params = ())
    for row in DBInterface.execute(db, String(sql), params)
        return row[1]
    end
    return nothing
end

function family_id(db::SQLite.DB, family_name::AbstractString)
    id = sqlite_first_value(db, "SELECT id FROM apn_family WHERE name = ?", (family_name,))
    id === nothing && error("Unknown APN family: $family_name")
    return id
end

function insert_apn_family_match!(db::SQLite.DB, function_id::Integer, match::APNFamilyMatch)
    name = apn_family_db_name(match.family)
    relation = match.exact ? "exact" : "equivalent"
    parameters = JSON3.write(match.parameters)
    DBInterface.execute(
        db,
        """
        INSERT OR IGNORE INTO apn_function_family (function_id, family_id, relation, parameters)
        VALUES (?, ?, ?, ?)
        """,
        (function_id, family_id(db, name), relation, parameters),
    )
end

function apn_function_id(db::SQLite.DB, dimension::Integer, local_id::Integer; required::Bool = false)
    id = sqlite_first_value(
        db,
        "SELECT id FROM apn_function WHERE dimension = ? AND local_id = ?",
        (dimension, local_id),
    )
    required && id === nothing && error("APN function dimension=$dimension local_id=$local_id does not exist")
    return id
end

function insert_apn_function_json!(db::SQLite.DB, json_source;
                                   id = nothing,
                                   source_label = nothing,
                                   notes = nothing,
                                   representation_type::AbstractString = "polynomial_json",
                                   ignore_existing::Bool = false)
    data = json_source isa AbstractString ? polynomial_json_data(json_source) : json_source
    canonical_json = canonical_polynomial_json(data)
    polynomial = build_polynomial_from_json(data)
    dimension = polynomial_json_dimension(data)
    display_expression = string(polynomial)
    normalized_support = normalized_support_from_json(data)

    local_id = id
    if local_id !== nothing && ignore_existing
        existing_id = apn_function_id(db, dimension, local_id)
        existing_id !== nothing && return Int(existing_id)
    end

    if local_id === nothing
        DBInterface.execute(
            db,
            """
            INSERT INTO apn_function
                (dimension, canonical_expression, canonical_json, normalized_support, source_label, notes)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (dimension, display_expression, canonical_json, normalized_support, source_label, notes),
        )
        function_id = SQLite.last_insert_rowid(db)
    else
        verb = ignore_existing ? "INSERT OR IGNORE" : "INSERT"
        DBInterface.execute(
            db,
            """
            $verb INTO apn_function
                (local_id, dimension, canonical_expression, canonical_json, normalized_support, source_label, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (local_id, dimension, display_expression, canonical_json, normalized_support, source_label, notes),
        )
        function_id = Int(apn_function_id(db, dimension, local_id; required = true))
    end

    DBInterface.execute(
        db,
        """
        INSERT OR IGNORE INTO apn_function_representation
            (function_id, representation_type, representation_json, display_expression, parameters)
        VALUES (?, ?, ?, ?, ?)
        """,
        (function_id, representation_type, canonical_json, display_expression, nothing),
    )

    for match in classify_family(symbolic_apn_function_from_json(data))
        insert_apn_family_match!(db, function_id, match)
    end

    return function_id
end

function apn_table_rows(table_json::AbstractString)
    header, rows = json_table(table_json)
    header == ["ID", "F(x)"] || error("APN table header must be [\"ID\", \"F(x)\"]")

    return rows
end

function insert_apn_function_table_json!(db::SQLite.DB, table_json::AbstractString;
                                         dimension::Integer,
                                         modulus::AbstractString,
                                         basis::AbstractString = "power",
                                         source_label = nothing,
                                         notes = nothing,
                                         ignore_existing::Bool = true)
    inserted_ids = Int[]

    SQLite.transaction(db) do
        for row in apn_table_rows(table_json)
            values = collect(row)
            length(values) == 2 || error("Each APN table row must have ID and F(x)")

            id = Int(values[1])
            expression = String(values[2])
            data = polynomial_expression_json(
                expression;
                dimension = dimension,
                modulus = modulus,
                basis = basis,
            )

            insert_apn_function_json!(
                db,
                data;
                id = id,
                source_label = source_label,
                notes = notes,
                representation_type = "table_expression_json",
                ignore_existing = ignore_existing,
            )
            push!(inserted_ids, id)
        end
    end

    return inserted_ids
end

function insert_apn_function_table_json!(path::AbstractString, table_json::AbstractString; kwargs...)
    return with_apn_database(path) do db
        insert_apn_function_table_json!(db, table_json; kwargs...)
    end
end

function insert_apn_functions!(db::SQLite.DB, json_sources;
                               source_label = nothing,
                               notes = nothing,
                               representation_type::AbstractString = "polynomial_json")
    inserted_ids = Int[]

    SQLite.transaction(db) do
        for json_source in json_sources
            function_id = insert_apn_function_json!(
                db,
                json_source;
                source_label = source_label,
                notes = notes,
                representation_type = representation_type,
            )
            push!(inserted_ids, Int(function_id))
        end
    end

    return inserted_ids
end

function insert_apn_functions!(path::AbstractString, json_sources; kwargs...)
    return with_apn_database(path) do db
        insert_apn_functions!(db, json_sources; kwargs...)
    end
end
