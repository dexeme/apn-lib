const APN_DB_DIR = normpath(@__DIR__)
const APN_PROJECT_ROOT = dirname(APN_DB_DIR)
const APN_JSON_DIR = joinpath(APN_DB_DIR, "json")
const APN_DB_PATH = joinpath(APN_DB_DIR, "apn_functions.sqlite")
const APN_SCHEMA_PATH = joinpath(APN_DB_DIR, "schema.sql")

apn_json_path(filename::AbstractString) = joinpath(APN_JSON_DIR, filename)
read_apn_json(filename::AbstractString) = read(apn_json_path(filename), String)
