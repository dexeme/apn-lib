using Test
using Random
using Nemo
using APNLib

const TEST_RNG = MersenneTwister(0xEA7)

identity_lut(n::Int)::Vector{Int} = collect(0:(space_size(n) - 1))

function gold_lut(n::Int)::Vector{Int}
    field = GF(2, n, "a")
    R, x = Nemo.polynomial_ring(field, "x")
    return univariate_to_lut(x^3, n)
end

function gf2_rank_matrix(matrix::Matrix{Int})::Int
    rows, cols = size(matrix)
    work = copy(matrix)
    rank = 0
    pivot_row = 1

    for col in 1:cols
        pivot = findfirst(row -> work[row, col] == 1, pivot_row:rows)
        pivot === nothing && continue

        pivot += pivot_row - 1
        work[pivot_row, :], work[pivot, :] = copy(work[pivot, :]), copy(work[pivot_row, :])

        for row in 1:rows
            if row != pivot_row && work[row, col] == 1
                work[row, :] .= xor.(work[row, :], work[pivot_row, :])
            end
        end

        rank += 1
        pivot_row += 1
        pivot_row > rows && break
    end

    return rank
end

function random_invertible_matrix(n::Int; rng::AbstractRNG = TEST_RNG)::Matrix{Int}
    while true
        matrix = rand(rng, 0:1, n, n)
        gf2_rank_matrix(matrix) == n && return matrix
    end
end

function matrix_apply(matrix::Matrix{Int}, x::Int, n::Int)::Int
    y = 0

    for row in 1:n
        bit = 0
        for col in 1:n
            bit = xor(bit, matrix[row, col] & ((x >> (col - 1)) & 1))
        end
        y |= bit << (row - 1)
    end

    return y
end

function linear_lut_from_matrix(matrix::Matrix{Int}, n::Int)::Vector{Int}
    return [matrix_apply(matrix, x, n) for x in 0:(space_size(n) - 1)]
end

function random_linear_permutation_lut(n::Int; rng::AbstractRNG = TEST_RNG)::Vector{Int}
    return linear_lut_from_matrix(random_invertible_matrix(n; rng = rng), n)
end

function cyclic_shift_lut(n::Int)::Vector{Int}
    matrix = zeros(Int, n, n)
    for col in 1:n
        matrix[mod1(col + 1, n), col] = 1
    end
    return linear_lut_from_matrix(matrix, n)
end

function compose_ea(A1::Vector{Int}, F::Vector{Int}, A2::Vector{Int}, A::Vector{Int}, n::Int)::Vector{Int}
    # Tested EA convention: G(x) = A1(F(A2(x))) + A(x) over GF(2)^n.
    return [xor(A1[F[A2[x + 1] + 1] + 1], A[x + 1]) for x in 0:(space_size(n) - 1)]
end

function validates_ea(eq::EAEquivalence, F::Vector{Int}, G::Vector{Int}, n::Int)::Bool
    return all(x -> G[x + 1] == xor(eq.L1[F[eq.A2[x + 1] + 1] + 1], eq.A[x + 1]),
               0:(space_size(n) - 1))
end

@testset "EA Equivalence Validation (Kaleyski Scalpel)" begin
    n = 7
    id = identity_lut(n)

    @testset "Test 1: Identity (Basic Sanity)" begin
        # The Gold function x^3 is APN and quadratic for odd n. When F = G, the
        # trivial equivalence is L1 = Id, A2 = Id, and A = 0. This sanity case
        # checks the LUT's 0-based mathematical indexing against Julia's 1-based
        # array indexing.
        F = gold_lut(n)
        G = copy(F)

        eq = first_ea_equivalence(F, G, n, is_quadratic = true)

        @test eq !== nothing
        @test eq.L1 == id
        @test eq.A2 == id
        @test affine_lut(eq.A, n)
        @test validates_ea(eq, F, G, n)
    end

    @testset "Test 2: Known Equivalence Reconstruction (The Disguise)" begin
        # Kaleyski/Heggebakk split the search into the external L1 map, the
        # internal A2 map, and the affine residual A. The implemented engine
        # normalizes F(0)=G(0)=0 and, in the quadratic path, tests A2 with zero
        # translation. We therefore use nontrivial linear maps, the affine
        # subcase with zero constant term.
        F = gold_lut(n)
        A1_secret = id
        A2_secret = cyclic_shift_lut(n)
        A_secret = id
        G = compose_ea(A1_secret, F, A2_secret, A_secret, n)

        @test affine_lut(A1_secret, n)
        @test affine_lut(A2_secret, n)
        @test affine_lut(A_secret, n)

        eq = first_ea_equivalence(F, G, n, is_quadratic = true)

        @test eq !== nothing
        @test affine_lut(eq.A, n)
        @test validates_ea(eq, F, G, n)
    end

    @testset "Test 3: Negative Case (Fast Partition Rejection)" begin
        # Algorithm 1 partitions values by the multiplicities of four-element
        # sums. These multiplicities are EA invariants. The linear identity has
        # a sum spectrum/signature incompatible with the Gold APN x^3 function,
        # so the external search must fail before the permutation DFS starts.
        F = gold_lut(n)
        G = identity_lut(n)

        partition_F = partition_by_multiplicity(F, n, k = 4)
        partition_G = partition_by_multiplicity(G, n, k = 4)
        signature_F = sort(zip(partition_F.multiplicities, length.(partition_F.blocks)) |> collect)
        signature_G = sort(zip(partition_G.multiplicities, length.(partition_G.blocks)) |> collect)

        @test signature_F != signature_G
        @test prepare_external_reconstruction(F, G, n, k = 4) === nothing
        @test first_ea_equivalence(F, G, n, is_quadratic = true) === nothing
    end
end
