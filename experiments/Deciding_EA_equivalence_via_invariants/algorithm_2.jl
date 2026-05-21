# Algorithm 2 Finding all linear permutations respecting a pair of partitions.
# --------------------------------------------------------------------------------
# Input : Two partitions K_i and C_i, a basis B, and a set U of possible images.
# Output: All linear permutations L1 such that L1(K_i) = C_i.
#
# Set L1(0) <- 0
# return assign(1)
#
# procedure assign(i)
#     if i = m + 1, return {L1}
#     Results <- empty
#     for c_i in U do
#         for x in Span({b1, ..., b_{i-1}}) do
#             L1(x + b_i) <- L1(x) + c_i
#             check that x + b_i and L1(x + b_i) are in corresponding blocks
#         end
#         if the partition is preserved, recurse
#     end
#     return Results

using APNLib

function kaleyski_assign_external_linear_maps!(results::Vector{Vector{Int}},
                                               K::KaleyskiMultiplicityPartition,
                                               C::KaleyskiMultiplicityPartition,
                                               B::Vector{Int},
                                               U::Vector{Int},
                                               n::Int,
                                               L1::Vector{Int},
                                               domain_span::Vector{Int},
                                               image_pivots::Vector{Int},
                                               i::Int)
    if i > n
        push!(results, copy(L1))
        return
    end

    basis_value = B[i]

    for c_i in U
        kaleyski_gf2_is_independent(c_i, image_pivots) || continue

        partition_preserved = true
        new_domain_values = Vector{Int}(undef, length(domain_span))
        new_image_values = Vector{Int}(undef, length(domain_span))

        for index in eachindex(domain_span)
            x = domain_span[index]
            y = xor(x, basis_value)
            y_image = xor(L1[x + 1], c_i)

            if C.block_index[y_image + 1] != K.block_index[y + 1]
                partition_preserved = false
                break
            end

            new_domain_values[index] = y
            new_image_values[index] = y_image
        end

        partition_preserved || continue

        old_span_length = length(domain_span)
        old_pivots = copy(image_pivots)
        kaleyski_gf2_add_pivot!(image_pivots, c_i)

        for index in eachindex(new_domain_values)
            L1[new_domain_values[index] + 1] = new_image_values[index]
            push!(domain_span, new_domain_values[index])
        end

        kaleyski_assign_external_linear_maps!(results, K, C, B, U, n, L1, domain_span, image_pivots, i + 1)

        resize!(domain_span, old_span_length)
        image_pivots .= old_pivots

        for value in new_domain_values
            L1[value + 1] = -1
        end
    end

    return
end

function algorithm_2(K::KaleyskiMultiplicityPartition,
                     C::KaleyskiMultiplicityPartition,
                     B::Vector{Int},
                     U::Vector{Int},
                     n::Int;
                     parallel::Bool = Threads.nthreads() > 1,
                     log_level = :quiet)
    APNLib.check_length(B, n, name = "basis")

    field_size = APNLib.space_size(n)
    L1 = fill(-1, field_size)
    L1[1] = 0

    results = Vector{Vector{Int}}()
    domain_span = [0]
    image_pivots = zeros(Int, n)

    kaleyski_assign_external_linear_maps!(results, K, C, B, U, n, L1, domain_span, image_pivots, 1)
    return results
end
