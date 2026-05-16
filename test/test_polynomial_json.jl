using Test
using DBInterface
using Nemo
using SQLite
using APNLib

@testset "Polynomial JSON representation" begin
    json_str = """
    {
      "field": "GF(2^7)",
      "basis": "power",
      "modulus": "t^7+t+1",
      "terms": [
        {"exp": 3, "coef": 22},
        {"exp": 9, "coef": 5}
      ]
    }
    """

    polynomial = build_polynomial_from_json(json_str)
    field = base_ring(parent(polynomial))

    @test coeff(polynomial, 3) == int_to_field_element(22, field, 7)
    @test coeff(polynomial, 9) == int_to_field_element(5, field, 7)
    @test coeff(polynomial, 4) == zero(field)

    expression_data = polynomial_expression_json(
        "g32x96 + gx33 + x3";
        dimension = 7,
        modulus = "t^7+t+1",
    )
    expression_polynomial = build_polynomial_from_json(expression_data)
    expression_field = base_ring(parent(expression_polynomial))

    @test coeff(expression_polynomial, 96) == gen(expression_field)^32
    @test coeff(expression_polynomial, 33) == gen(expression_field)
    @test coeff(expression_polynomial, 3) == one(expression_field)
end

@testset "APN function batch insert" begin
    json_str = """
    {
      "field": "GF(2^7)",
      "basis": "power",
      "modulus": "t^7+t+1",
      "terms": [
        {"exp": 3, "coef": 1},
        {"exp": 9, "coef": 1},
        {"exp": 18, "coef": 1},
        {"exp": 36, "coef": 1},
        {"exp": 72, "coef": 1},
        {"exp": 17, "coef": 1},
        {"exp": 34, "coef": 1},
        {"exp": 68, "coef": 1}
      ]
    }
    """

    db_path = tempname() * ".sqlite"
    ids = insert_apn_functions!(db_path, [json_str, json_str]; source_label = "test")
    @test length(ids) == 2

    db = SQLite.DB(db_path)
    try
        function_count = first(DBInterface.execute(db, "SELECT COUNT(*) AS count FROM apn_function"))[1]
        representation_count = first(DBInterface.execute(db, "SELECT COUNT(*) AS count FROM apn_function_representation"))[1]
        family_count = first(DBInterface.execute(db, "SELECT COUNT(*) AS count FROM apn_function_family"))[1]

        @test function_count == 2
        @test representation_count == 2
        @test family_count >= 2
    finally
        SQLite.close(db)
        isfile(db_path) && rm(db_path)
    end
end

@testset "APN function table JSON insert" begin
    table_json = """
    [["ID","F(x)"],
     [1,"x3"],
     [19,"g32x96 + gx33 + x3"]]
    """

    db_path = tempname() * ".sqlite"
    ids = insert_apn_function_table_json!(
        db_path,
        table_json;
        dimension = 7,
        modulus = "t^7+t+1",
        source_label = "table-json-test",
    )

    @test ids == [1, 19]

    db = SQLite.DB(db_path)
    try
        function_count = first(DBInterface.execute(db, "SELECT COUNT(*) AS count FROM apn_function"))[1]
        representation_json = first(DBInterface.execute(
            db,
            """
            SELECT r.representation_json
            FROM apn_function_representation r
            JOIN apn_function f ON f.id = r.function_id
            WHERE f.dimension = ? AND f.local_id = ?
            """,
            (7, 19),
        ))[1]

        @test function_count == 2
        @test occursin("\"coef_power\":32", representation_json)
        @test occursin("\"coef_power\":1", representation_json)
    finally
        SQLite.close(db)
        isfile(db_path) && rm(db_path)
    end
end

@testset "APN invariant table insert" begin
    functions_json = """
    [["ID","F(x)"],
     [1,"x3"],
     [2,"x9"]]
    """
    gamma_json = """
    [["Gamma-rank","ID"],
     [3610,"1,2"]]
    """
    delta_json = """
    [["Delta-rank","ID"],
     [198,1],
     [210,2]]
    """
    multiplier_json = """
    [["Multiplier group order","ID"],
     [113792,1],
     [128,2]]
    """

    db_path = tempname() * ".sqlite"
    insert_apn_function_table_json!(
        db_path,
        functions_json;
        dimension = 7,
        modulus = "t^7+t+1",
    )

    invariant_columns = Dict(
        "Gamma-rank" => "gamma_rank",
        "Delta-rank" => "delta_rank",
        "Multiplier group order" => "multiplier_group_order",
    )
    updated_ids = with_apn_database(db_path) do db
        ids = Int[]
        SQLite.transaction(db) do
            for table_json in [gamma_json, delta_json, multiplier_json]
                header, rows = json_table(table_json)
                require_json_columns(header, ["ID"])
                value_column = only(setdiff(header, ["ID"]))
                target_column = invariant_columns[value_column]

                for row in rows
                    data = json_row_dict(header, row)
                    for local_id in split_int_values(data["ID"])
                        function_id = apn_function_id(db, 7, local_id; required = true)
                        upsert_table_values!(
                            db,
                            "apn_invariant",
                            "function_id",
                            apn_public_function_id(7, local_id),
                            Dict("dimension" => 7, target_column => Int(data[value_column])),
                        )
                        push!(ids, local_id)
                    end
                end
            end
        end
        ids
    end

    @test length(updated_ids) == 6

    db = SQLite.DB(db_path)
    try
        row = first(DBInterface.execute(
            db,
            """
            SELECT dimension, gamma_rank, delta_rank, multiplier_group_order
            FROM apn_invariant
            WHERE function_id = ?
            """,
            ("7.1",),
        ))

        @test row[1] == 7
        @test row[2] == 3610
        @test row[3] == 198
        @test row[4] == 113792
    finally
        SQLite.close(db)
        isfile(db_path) && rm(db_path)
    end
end

@testset "APN equivalence and Walsh spectrum insert" begin
    functions_json = """
    [["ID","F(x)"],
     [1,"x3"],
     [2,"x9"]]
    """
    equivalence_json = """
    [["N∘","Equivalent to","Walsh spectrum"],
     [7.1,"Gold","Gold"],
     [7.2,"C3","Gold"]]
    """

    db_path = tempname() * ".sqlite"
    insert_apn_function_table_json!(
        db_path,
        functions_json;
        dimension = 7,
        modulus = "t^7+t+1",
    )

    updated_ids = with_apn_database(db_path) do db
        ids = Int[]
        header, rows = json_table(equivalence_json)
        require_json_columns(header, ["N∘", "Equivalent to", "Walsh spectrum"])

        SQLite.transaction(db) do
            for row in rows
                data = json_row_dict(header, row)
                parts = split(string(data["N∘"]), ".")
                dimension = parse(Int, parts[1])
                local_id = parse(Int, parts[2])
                function_id = apn_function_id(db, dimension, local_id; required = true)
                update_table_values!(
                    db,
                    "apn_function",
                    Dict("equivalent_to" => String(data["Equivalent to"]));
                    where_column = "id",
                    where_value = function_id,
                )
                upsert_table_values!(
                    db,
                    "apn_invariant",
                    "function_id",
                    apn_public_function_id(dimension, local_id),
                    Dict(
                        "dimension" => dimension,
                        "walsh_spectrum" => String(data["Walsh spectrum"]),
                    ),
                )
                push!(ids, Int(function_id))
            end
        end
        ids
    end
    @test updated_ids == [1, 2]

    db = SQLite.DB(db_path)
    try
        first_row = first(DBInterface.execute(
            db,
            """
            SELECT f.equivalent_to, i.walsh_spectrum
            FROM apn_function f
            JOIN apn_invariant i ON i.function_id = CAST(f.dimension AS TEXT) || '.' || CAST(f.local_id AS TEXT)
            WHERE f.dimension = ? AND f.local_id = ?
            """,
            (7, 1),
        ))
        tenth_row = first(DBInterface.execute(
            db,
            """
            SELECT f.equivalent_to, i.walsh_spectrum
            FROM apn_function f
            JOIN apn_invariant i ON i.function_id = CAST(f.dimension AS TEXT) || '.' || CAST(f.local_id AS TEXT)
            WHERE f.dimension = ? AND f.local_id = ?
            """,
            (7, 2),
        ))

        @test first_row[1] == "Gold"
        @test first_row[2] == "Gold"
        @test tenth_row[1] == "C3"
        @test tenth_row[2] == "Gold"
    finally
        SQLite.close(db)
        isfile(db_path) && rm(db_path)
    end
end
