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

@testset "APN Search Optional CSV Output" begin
    identity1 = reshape([1], 1, 1)
    results_filename = tempname()

    solutions = APNSearch(
        1,
        identity1,
        identity1,
        max_solutions = 1,
        on_solution = sbox -> nothing,
        save_results = true,
        results_filename = results_filename,
    )

    lines = readlines(results_filename)

    @test solutions == [[0, 1]]
    @test lines[1] == "s0,s1,polynomial"
    @test lines[2] == "0,1,\"x^1\""
end
