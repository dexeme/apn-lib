using Test
using APNLib

function _gf8_mul(a::Int, b::Int)::Int
    result = 0
    aa = a
    bb = b

    for _ in 1:3
        if isodd(bb)
            result ⊻= aa
        end

        bb >>= 1
        carry = (aa & 0x04) != 0
        aa = (aa << 1) & 0x07

        if carry
            aa ⊻= 0x03
        end
    end

    return result
end

_gf8_square(a::Int)::Int = _gf8_mul(a, a)
_gf8_cube(a::Int)::Int = _gf8_mul(_gf8_square(a), a)

@testset "Orthoderivative" begin
    lut = [_gf8_cube(x) for x in 0:7]
    pi = orthoderivative(lut)

    @test pi[1] == 0
    @test pi == collect(0:7)

    for alpha in 1:7
        y = pi[alpha + 1]

        @test 1 <= y <= 7
        @test lut[(alpha ⊻ y) + 1] ⊻ lut[alpha + 1] ⊻ lut[y + 1] ⊻ lut[1] == 0
    end

    x9_lut = [_gf8_square(x) for x in 0:7]
    @test orthoderivative(x9_lut) == [0, 1, 1, 1, 1, 1, 1, 1]
end
