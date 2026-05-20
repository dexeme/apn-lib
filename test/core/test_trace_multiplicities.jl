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
    @test walsh_spectrum(identity_lut, 2) == vec(walsh_table)
    @test extended_walsh_spectrum(identity_lut, 2) == sort(abs.(vec(walsh_table)))
end

@testset "Fast Walsh Table Matches Direct Formula" begin
    n = 3
    field = GF(2, n, "g")
    gold_polynomial = APNFunction(monomial_expr(3))
    lut = univariate_to_lut(gold_polynomial, n)
    inputs = APNLib.field_elements(field, n)
    values = APNLib.function_values_to_field(lut, field, n)
    walsh_table = walsh_coefficient_table(lut, n)

    for a in 0:(space_size(n) - 1)
        for b in 0:(space_size(n) - 1)
            direct = APNLib.walsh_coefficient(values, inputs, inputs[a + 1], inputs[b + 1])
            @test walsh_table[a + 1, b + 1] == direct
        end
    end
end

@testset "Walsh Multiplicities" begin
    identity_lut = [0, 1, 2, 3]

    @test multiplicities_sigma(identity_lut, 2, 2) == Dict(0 => 4, 1 => 0, 2 => 0, 3 => 0)

    n = 3
    field = GF(2, n, "g")
    gold_polynomial = APNFunction(monomial_expr(3))
    lut = univariate_to_lut(gold_polynomial, n)
    elements = APNLib.field_elements(field, n)
    walsh_table = walsh_coefficient_table(lut, n)
    divisor = space_size(n)^2
    direct_multiplicities = Dict{Int, Int}()

    for (s_index, s) in pairs(elements)
        total = 0

        for a_index in eachindex(elements)
            for (b_index, b) in pairs(elements)
                total += trace_sign(b * s) * walsh_table[a_index, b_index]^4
            end
        end

        direct_multiplicities[s_index - 1] = div(total, divisor)
    end

    @test multiplicities_sigma(lut, n, 4) == direct_multiplicities
    @test multiplicities_sigma(gold_polynomial, n, 4) == direct_multiplicities

    dimensioned_gold_polynomial = APNFunction(n, monomial_expr(3))

    @test multiplicities_sigma(dimensioned_gold_polynomial, n, 4) == direct_multiplicities
end
