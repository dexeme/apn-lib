using APNLib
using SQLite

include(joinpath(@__DIR__, "..", "config.jl"))

const INVARIANT_COLUMN_MAP = Dict(
    "Gamma-rank" => "gamma_rank",
    "Delta-rank" => "delta_rank",
    "Multiplier group order" => "multiplier_group_order",
)

const INVARIANT_JSON_FILES = [
    "apn_gamma_rank_7.json",
    "apn_delta_rank_7.json",
    "apn_multiplier_group_order_7.json",
]

function insert_invariant_json!(db, table_json::AbstractString; dimension::Integer, id_column::AbstractString = "ID")
    header, rows = json_table(table_json)
    length(header) == 2 || error("Invariant JSON must have two columns")
    require_json_columns(header, [id_column])

    value_column = only(setdiff(header, [id_column]))
    target_column = get(INVARIANT_COLUMN_MAP, value_column, nothing)
    target_column !== nothing || error("Unsupported invariant column: $value_column")

    updated_ids = Int[]
    for row in rows
        data = json_row_dict(header, row)
        value = Int(data[value_column])

        for local_id in split_int_values(data[id_column])
            function_id = apn_function_id(db, dimension, local_id; required = true)
            upsert_table_values!(
                db,
                "apn_invariant",
                "function_id",
                function_id,
                Dict(target_column => value),
            )
            push!(updated_ids, local_id)
        end
    end

    return updated_ids
end

function insert_invariant_jsons!(db, table_jsons; dimension::Integer)
    updated_ids = Int[]
    SQLite.transaction(db) do
        for table_json in table_jsons
            append!(updated_ids, insert_invariant_json!(db, table_json; dimension = dimension))
        end
    end
    return updated_ids
end

updated_ids = with_apn_database(APN_DB_PATH) do db
    insert_invariant_jsons!(db, read_apn_json.(INVARIANT_JSON_FILES); dimension = 7)
end

println("Updated $(length(updated_ids)) invariant values for APN functions")
println("Database: $(APN_DB_PATH)")
