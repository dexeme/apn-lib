#Algorithm 1: Reconstructing the outer permutation L1 in L1 ◦ F ◦ L2 = G
#-----------------------------------------------------------------------
#Input : Two (n,m)-functions F and G
#Output: All linear permutations L1 of (F_2)^m respecting the partitions induced by F and G
#
#Partition (F_2)^m = B^F_0 ∪ · · · ∪ B^F_{2^m} = B^G_0 ∪ · · · ∪ B^G_{2^m} ;
#if (∃i)(#B^F_i != #B^G_i) then return ∅ end
#
#Let B = {b1, b2, . . . , bm} be a basis of (F_2)^m ;
#return Guess(B, 1, ∅) # recursively guess the values of L1 on B
#
#Function Guess(B, i, L):
#if i = n + 1 then
#    reconstruct L1 from its values on B ;
#    return L ∪ {L1} ;
#end
#
#Let j be such that b_i ∈ B^F_j ;
#for y ∈ B^G_j do
#    L1(b_i) ← y ;
#    # Check all currently known values of L1 for contradiction
#    contradiction ← false ;
#    for x ∈ Span(b1, b2, . . . , b_i) do
#        let j, k be such that x ∈ B^F_j, L1(x) ∈ B^G_k ;
#        if k != j then
#            contradiction ← true ;
#            break ;
#        end
#    end
#    if contradiction = false then
#        L ← L ∪ Guess(B, i + 1, L) ;
#    end
#end
#return L ;

