using Test
using Nemo
using APNLib

@testset "Dynamic DDT Updates" begin
    sbox = fill(-1, 4)
    ddt = zeros(Int, 4, 4)

    sbox[0 + 1] = 0
    @test addDDTInformation(0, sbox, ddt) == true
    @test all(ddt .== 0)

    sbox[3 + 1] = 3
    @test addDDTInformation(3, sbox, ddt) == true
    @test ddt[3 + 1, 3 + 1] == 2
    @test removeDDTInformation(3, sbox, ddt) == true
    @test all(ddt .== 0)
end

@testset "S-box Search State Helpers" begin
    sbox = [0, -1, -1, 3]

    @test isComplete(sbox) == false
    @test nextFreePosition(sbox) == 1
    @test isComplete([0, 1, 2, 3]) == true
end

@testset "S-box Polynomial Interpolation" begin
    identity_lut = collect(0:7)

    @test string(int_to_field_element(5, Nemo.GF(2, 3, "g"), 3)) == "g^2 + 1"
    @test format_sbox_polynomial(identity_lut, 3) == "x^1"
end

@testset "APN Search Matrix Constant Output" begin
    identity3 = Int[1 0 0; 0 1 0; 0 0 1]
    constants_filename = joinpath("tuples", "AllTuplesMatrices3.jl")
    original_text = read(constants_filename, String)

    try
        solutions = APNSearch(
            3,
            identity3,
            identity3,
            max_solutions = 1,
            on_solution = sbox -> nothing,
            save_results = true,
            class_index = 1,
        )

        text = read(constants_filename, String)
        expected_line = "const ALL_TUPLES_3_1_SEARCH = Int[$(join(solutions[1], ", "))]"

        @test length(solutions) == 1
        @test occursin(expected_line, text)
    finally
        write(constants_filename, original_text)
    end
end
