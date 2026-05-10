using APNLib
include(joinpath(@__DIR__, "..", "config.jl"))

json_path = apn_json_path("apn_dimension_7.json")

table_json = read_apn_json("apn_dimension_7.json")

inserted_ids = insert_apn_function_table_json!(
    APN_DB_PATH,
    table_json;
    dimension = 7,
    modulus = "t^7+t+1",
    source_label = "apn_dimension_7",
)

println("Inserted or reused $(length(inserted_ids)) APN functions from $(json_path)")
println("Database: $(APN_DB_PATH)")
