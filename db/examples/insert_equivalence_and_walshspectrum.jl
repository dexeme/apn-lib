using APNLib
using SQLite

include(joinpath(@__DIR__, "..", "config.jl"))

const EQUIVALENCE_JSON_FILE = "apn_equivalence_and_walshspectrum.json"
const FUNCTION_COLUMN_MAP = Dict(
    "Equivalent to" => "equivalent_to",
    "Walsh spectrum" => "walsh_spectrum",
)

function dimension_local_id(value, last_local_ids::Dict{Int, Int})
    parts = split(string(value), ".")
    length(parts) == 2 || error("Expected ID in dimension.local_id form, got $value")

    dimension = parse(Int, parts[1])
    parsed_local_id = parse(Int, parts[2])
    previous_local_id = get(last_local_ids, dimension, 0)
    local_id = parsed_local_id <= previous_local_id ? previous_local_id + 1 : parsed_local_id
    last_local_ids[dimension] = local_id
    return dimension, local_id
end

function insert_function_column_json!(db, table_json::AbstractString;
                                      id_column::AbstractString = "ID",
                                      column_map = FUNCTION_COLUMN_MAP,
                                      skip_missing_function::Bool = true)
    header, rows = json_table(table_json)
    require_json_columns(header, [id_column; collect(keys(column_map))])

    updated_ids = Int[]
    last_local_ids = Dict{Int, Int}()

    SQLite.transaction(db) do
        for row in rows
            data = json_row_dict(header, row)
            dimension, local_id = dimension_local_id(data[id_column], last_local_ids)
            function_id = apn_function_id(db, dimension, local_id)

            if function_id === nothing
                skip_missing_function && continue
                error("APN function dimension=$dimension local_id=$local_id does not exist")
            end

            values = Dict(target_column => String(data[source_column])
                          for (source_column, target_column) in column_map)
            update_table_values!(
                db,
                "apn_function",
                values;
                where_column = "id",
                where_value = function_id,
            )
            push!(updated_ids, Int(function_id))
        end
    end

    return updated_ids
end

json_path = apn_json_path(EQUIVALENCE_JSON_FILE)
updated_ids = with_apn_database(APN_DB_PATH) do db
    insert_function_column_json!(db, read_apn_json(EQUIVALENCE_JSON_FILE))
end

println("Updated equivalence and Walsh spectrum for $(length(updated_ids)) APN functions from $(json_path)")
println("Database: $(APN_DB_PATH)")
