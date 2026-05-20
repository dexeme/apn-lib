#Algorithm 2 Finding all linear permutations respecting a pair of partitions.
#-- -- -- -- -- -- -- -- ---- -- -- -- -- -- -- -- -- ---- -- -- -- -- -- -- -- -- ---- -- -- -- -- -- -- -- --
#Input:Two partitions (F_2)^m = K_1 ⊕ K_2 ⊕ · · · ⊕ K_s and (F_2)^m = C_1 ⊕ C_2 ⊕ · · · ⊕ C_s of the vector space (F_2)^m, a basis B = {b1, b2, · · ·, bm} of (F_2)^m, and a set U of possible values for the images of B
#Output:All linear permutations L1 of (F_2)^m such that L1(K_i) = C_i for 1 ≤ i ≤ s
#
#Set L1(0) ← 0;
#return assign(1)
#
#procedure assign(i) ;
#if i = m + 1 then
#
#return {L1} ;
#end
#Results ← ∅ ;
#for c_i ∈ U do
#    partitionPreserved ← true ;
#    for x ∈ Span({b1, · · · , b_{i−1}}) do
#        L1(x + b_i) ← L1(x) + c_i ;
#        find j such that x + b_i ∈ K_j ;
#        if L1(x + b_i) ∉ C_j then
#            partitionPreserved ← false ;
#            break ;
#        end
#    end
#    if partitionPreserved then
#        Results ← Results ∪ assign(i + 1) ;
#    end
#end
#return Results# algorithm2.jl
## Julia Script
#
##=
#Description:
#Author: tiago
#Date: 20/05/2026
#=#
#
#function main()
#    println("Hello, Julia!")
#end
#
#main()
