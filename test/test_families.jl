using Test
using APNLib

@testset "APN family recognition" begin
    c4 = APNFunction(x(3), Tr(7, x(9)))
    c3 = APNFunction(8, x(3), x(17), pterm(48, 18), pterm(3, 33), pterm(34), x(48))
    c5 = APNFunction(x(3), Tr(39, x(9), x(18)))

    @test belongs_to_family(c4, :C4)
    @test belongs_to_family(c3, :C3)
    @test belongs_to_family(c5, :C5)

    c3_match = only(filter(match -> match.family == :C3, classify_family(c3)))
    @test c3_match.family == :C3
    @test c3_match.parameters[:m] == 4
    @test c3_match.parameters[:q] == 16
    @test c3_match.parameters[:i] == 1
    @test c3_match.exact == false
end
