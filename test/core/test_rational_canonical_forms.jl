using Test
using APNLib

@testset "RCF Block Construction" begin
    F = gf2()
    blocks = blocks_for_rcf(3)

    @test length(blocks) == 4
    @test all(nrows(block) == 3 && ncols(block) == 3 for block in blocks)
    @test all(base_ring(block) == F for block in blocks)
    @test all(!iszero(det(block)) for block in blocks)
    @test_throws ErrorException blocks_for_rcf(0)
    @test_throws ErrorException APNLib.companion_matrix_gf2(F, [F(0), F(1)])
end

@testset "RCF Enumeration Counts" begin
    @test length(get_rcfs(2)) == 3
    @test length(get_rcfs(3)) == 6
    @test length(get_rcfs(4)) == 14
    @test_throws ErrorException get_rcfs(0)
end

@testset "Matrix Similarity Invariants" begin
    F = gf2()
    identity2 = identity_matrix(F, 2)
    jordan2 = matrix(F, 2, 2, [1, 1, 0, 1])
    swap = matrix(F, 2, 2, [0, 1, 1, 0])

    @test matrix_is_similar(identity2, identity2)
    @test matrix_is_similar(jordan2, swap)
    @test !matrix_is_similar(identity2, swap)
    @test_throws ErrorException matrix_is_similar(identity2, identity_matrix(F, 3))
end
