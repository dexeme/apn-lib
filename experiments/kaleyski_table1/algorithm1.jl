#Algorithm 1 General framework for reconstructing the outer permutation.
#-----------------------------------------------------------------------
#Input : Two (n, m)-functions F and G
#Output: All linear permutations L1 of (F_2)^m respecting the partitions induced by F and G
#
#for s ∈ (F_2)^m do
#    compute the number M^F_4(0, s) of (x1, x2, x3, x1 + x2 + x3) ∈ T_4(0) such that
#        F(x1) + F(x2) + F(x3) + F(x1 + x2 + x3) = s ;
#    compute the number M^G_4(0, s) of (x1, x2, x3, x1 + x2 + x3) ∈ T_4(0) such that
#        G(x1) + G(x2) + G(x3) + G(x1 + x2 + x3) = s ;
#end
#
#partition (F_2)^m = K_1 ⊕ K_2 ⊕ · · · ⊕ K_s so that M^F_4(0, s1) = M^F_4(0, s2) for s1 ∈ K_i and s2 ∈ K_j if and only if i = j ;
#partition (F_2)^m = C_1 ⊕ C_2 ⊕ · · · ⊕ C_s' so that M^G_4(0, s1) = M^G_4(0, s2) for s1 ∈ C_i and s2 ∈ C_j if and only if i = j ;
#
#if s != s' then
#    return ∅
#end
#
#rearrange C_1, C_2, . . . , C_s if necessary so that M^F_4(0, s1) = M^G_4(0, s2) where s1 ∈ K_i, s2 ∈ C_i for any 1 ≤ i ≤ s ;
#
#if #C_i != #K_i for some i in 1 ≤ i ≤ s then
#    return ∅
#end