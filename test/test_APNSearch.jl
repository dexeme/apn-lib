using Test
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
