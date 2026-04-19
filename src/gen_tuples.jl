# this is to generate all the LE-automorphisms to consider for n = 6,7,8, listed in Corollary 1
# Dependencies:
#import Pkg;
#Pkg.add("Combinatorics")

using LinearAlgebra
using Combinatorics
include("check_order_space.jl")
include("matrix_power_gf2.jl")


# converts the matrix A to a look-up table
function matrix_to_sbox(A::Matrix)
    T = []
    nrows = size(A, 1)
    # Generate all vectors in GF(2)^nrows
    for v_int in 0:(2^nrows - 1)
        v = digits(v_int, base=2, pad=nrows)
        result = (A * v) .% 2
        # Convert binary vector to integer
        result_int = sum(result[i] * 2^(i-1) for i in 1:length(result))
        push!(T, result_int)
    end
    return T
end

# checks if a matrix tuple (A,B) can be a self invariance for a permutation
# i.e., check whether, for all i, Ord(A,i) has the same dimension as Ord(B,i)
function is_permutation_tuple(A::Matrix, B::Matrix)
    for i in 0:(2^(size(A, 1)) - 1)
        if check_order_space(A, i) != check_order_space(B, i)
            return false
        end
    end
    return true
end

# returns all invertible companion matrices of dimension d over GF(2)
function blocks_for_rcf(d::Int)
    R = []
    # Generate all binary vectors of length d+1 where first element is 1
    for v_int in 0:(2^(d+1) - 1)
        v = digits(v_int, base=2, pad=d+1)
        if v[end] == 1 && v[1] == 1  # Check first and last elements
            # Create companion matrix from polynomial coefficients
            C = companion_matrix_gf2(v)
            push!(R, C)
        end
    end
    return R
end

# returns all rational canonical forms of GF(2)-matrices of dimension n
function get_rcfs(n::Int)
    R = []
    for partition in partitions(n)
        V = []
        for block_size in partition
            push!(V, blocks_for_rcf(block_size))
        end
        # Cartesian product of all block combinations
        for combo in Iterators.product(V...)
            app = true
            poly = minimal_polynomial_gf2(combo[1])
            for i in 1:length(combo)
                if !divides(minimal_polynomial_gf2(combo[i]), poly)
                    app = false
                end
                poly = minimal_polynomial_gf2(combo[i])
            end
            if app
                # Create block diagonal matrix
                BDM = block_diagonal_matrix_gf2(reverse(collect(combo)))
                push!(R, BDM)
            end
        end
    end
    return R
end

# checks if A=[A[0],A[1]] is power similar to B=[B[0],B[1]]
function is_power_similar(A::Vector, B::Vector)
    ord_A0 = multiplicative_order_gf2(A[1])
    ord_A1 = multiplicative_order_gf2(A[2])
    ord_B0 = multiplicative_order_gf2(B[1])
    ord_B1 = multiplicative_order_gf2(B[2])

    if ord_A0 == ord_B0 && ord_A1 == ord_B1
        for i in 0:(max(ord_A0, ord_A1) - 1)
            A0_power = matrix_power_gf2(A[1], i)
            A1_power = matrix_power_gf2(A[2], i)
            if is_similar_gf2(A0_power, B[1]) && is_similar_gf2(A1_power, B[2])
                return true
            end
        end
    end
    return false
end

# checks if A=[A[0],A[1]] is extended power similar to B=[B[0],B[1]]
function is_extended_power_similar(A::Vector, B::Vector)
    if is_power_similar(A, B)
        return true
    end
    if is_power_similar([inv_gf2(A[1]), inv_gf2(A[2])], [B[2], B[1]])
        return true
    end
    return false
end

# generates all possible matrix tuples that need to be considered for self equivalence
# (i.e., 17 for n=6, 27 for n=7, and 32 for n=8)
function gen_permutation_tuples(n::Int)
    T = get_rcfs(n)
    RC = []
    for t in T
        if is_prime(multiplicative_order_gf2(t))
            push!(RC, t)
        end
    end

    T = collect(Iterators.product(RC, RC))
    G = []
    for t in T
        if multiplicative_order_gf2(t[1]) == multiplicative_order_gf2(t[2])
            push!(G, [t[1], t[2], 1])
        end
    end

    GG = []
    ctr = 0
    for i in 1:length(G)
        for j in 1:length(G)
            if i != j
                if G[j][3] == 1
                    if is_extended_power_similar(G[j], G[i])
                        G[i][3] = 0
                        ctr += 1
                        break
                    end
                end
            end
        end
        println("$i, $ctr")
    end

    for g in G
        if g[3] == 1
            push!(GG, g[1:2])
        end
    end

    GGG = []
    for i in 1:length(GG)
        t = GG[i]
        println("$i, $(length(GGG))")
        if is_permutation_tuple(t[1], t[2])
            push!(GGG, t)
        end
    end

    return GGG
end

# generates the file AllTuples.h for the matrix tuples in list T
function generate_tuples_file(T::Vector)
    L = []
    for i in 1:length(T)
        append!(L, [vcat(matrix_to_sbox(T[i][1]), matrix_to_sbox(T[i][2]))])
    end

    open("AllTuples.h", "w") do file
        write(file, "#define N_TUPLES $(length(L))\n\n")
        write(file, "int AllTuples[$(length(L))][$(length(L[1]))]={\n")

        for ind in 1:(length(L) - 1)
            write(file, "{")
            for i in 1:(length(L[ind]) - 1)
                write(file, "$(L[ind][i]), ")
            end
            write(file, "$(L[ind][end])},\n")
        end

        write(file, "{")
        for i in 1:(length(L[end]) - 1)
            write(file, "$(L[end][i]), ")
        end
        write(file, "$(L[end][end])}\n")
        write(file, "};\n")
    end
end

# Main execution
T6 = gen_permutation_tuples(6)
T7 = gen_permutation_tuples(7)
T8 = gen_permutation_tuples(8)