using Test
using Nemo
using APNLib

@testset "Absolute Trace and Walsh Spectrum" begin
    field = GF(2, 2, "g")
    elements = field_elements(field, 2)

    @test absolute_trace_bit(zero(field)) == 0
    @test all(absolute_trace_bit(element) in (0, 1) for element in elements)

    identity_lut = [0, 1, 2, 3]
    walsh_table = walsh_coefficient_table(identity_lut, 2)

    @test all(walsh_table[a, b] == (a == b ? 4 : 0) for a in 1:4, b in 1:4)
    @test walsh_spectrum(identity_lut, 2) == vec(walsh_table')
end

@testset "Walsh Multiplicities" begin
    identity_lut = [0, 1, 2, 3]

    @test multiplicities_sigma(identity_lut, 2, 2) == Dict(0 => 4, 1 => 0, 2 => 0, 3 => 0)
end

