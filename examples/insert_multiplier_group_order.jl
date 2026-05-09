using APNLib

project_root = normpath(joinpath(@__DIR__, ".."))
json_path = joinpath(project_root, "json", "apn_multiplier_group_order_7.json")
db_path = joinpath(project_root, "db", "apn_functions.sqlite")

table_json = read(json_path, String)

updated_ids = insert_invariants(
    db_path,
    [table_json];
    dimension = 7,
)

println("Updated Multiplier Group Order for $(length(updated_ids)) APN functions from $(json_path)")
println("Database: $(db_path)")
