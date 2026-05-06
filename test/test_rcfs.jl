using Test
using APNLib

@testset "RCF Block Construction" begin
    blocks = blocks_for_rcf(3)

    @test length(blocks) == 4
end

@testset "RCF Enumeration Counts" begin
    @test length(get_rcfs(2)) == 3
    @test length(get_rcfs(3)) == 6
    @test length(get_rcfs(4)) == 14
end
