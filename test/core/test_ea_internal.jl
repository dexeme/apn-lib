using Test
using Nemo
using APNLib

function field_multiplier_lut(multiplier, n::Int)::Vector{Int}
    field = parent(multiplier)
    lookup = APNLib.field_element_lookup(field, n)
    lut = Vector{Int}(undef, space_size(n))

    for x in 0:(space_size(n) - 1)
        y = multiplier * int_to_field_element(x, field, n)
        lut[x + 1] = lookup[y]
    end

    return lut
end

@testset "EA Internal Affine Reconstruction Algorithm 3" begin
    n = 6
    field = GF(2, n, "a")
    R, z = Nemo.polynomial_ring(field, "x")
    alpha = gen(field)

    F = univariate_to_lut(z^3, n)
    L2_secret = field_multiplier_lut(alpha^5, n)
    A_secret = field_multiplier_lut(alpha^9, n)
    G = [xor(F[L2_secret[x + 1] + 1], A_secret[x + 1]) for x in 0:(space_size(n) - 1)]

    @test F[1] == 0
    @test G[1] == 0

    solutions = reconstruct_internal_affine_maps(F, G, n, is_quadratic = true)

    @test (L2_secret, A_secret) in solutions
end
