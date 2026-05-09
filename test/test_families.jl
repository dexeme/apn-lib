using Test
using APNLib

@testset "APN family recognition" begin
    @test belongs_to_family("x3+Tr7(x9)", :C4)
    @test belongs_to_family("x3+x17+p48x18+p3x33+px34+x48", :C3; n = 8)
    @test belongs_to_family("x3+Tr39(x9+x18)", :C5)

    c3_match = only(filter(match -> match.family == :C3, classify_family("x3+x17+p48x18+p3x33+px34+x48"; n = 8)))
    @test c3_match.family == :C3
    @test c3_match.parameters[:m] == 4
    @test c3_match.parameters[:q] == 16
    @test c3_match.parameters[:i] == 1
    @test c3_match.exact == false
end
