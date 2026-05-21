#Algorithm 3 Reconstructing the inner permutation A2.
#----------------------------------------------------
#Input : Two (n, m)-functions F and G with F(0) = G(0) = 0
#Output: All affine permutations A2 of (F_2)^m such that F ◦ A2 + G is affine
#
#for x ∈ (F_2)^n do
#    D(x) ← (F_2)^n (initialize domains) ;
#    O^F_3(x) ← ∅ ;
#end
#for (x1, x2) ∈ ((F_2)^m)^2 do
#    t ← F(x1) + F(x2) + F(x1 + x2) ;
#    O^F_3(t) ← O^F_3(t) ∪ (x1, x2, x1 + x2) ;
#end
#for (x1, x2) ∈ ((F_2)^m)^2 do
#    t ← G(x1) + G(x2) + G(x1 + x2) ;
#    D(x1) ← D(x1) ∩ O^F_3(t) ;
#    D(x2) ← D(x2) ∩ O^F_3(t) ;
#    D(x1 + x2) ← D(x1 + x2) ∩ O^F_3(t) ;
#end
#Order the elements x ∈ (F_2)^m into x_i for 1 ≤ i ≤ 2^n, so that i < j =⇒ #D(x_i) ≤ #D(x_j) ;
#B ← ∅ (basis) ;
#Results ← ∅ ;
#for i = 1, 2, . . . , 2^m do
#    if x_i ∉ Span(B) then
#        B ← B ∪ {x_i} ;
#        if #B = m then
#            break ;
#        end
#    end
#end
#for c2 ∈ (F_2)^m do
#    for (v1, v2, . . . , vn) ∈ ∏_{x∈B} D(x) do
#        Let L2 be linear with L2(b_i) = v_i for B = (b1, . . . , bm) ;
#        A2 ← L2 + c2 ;
#        A ← F ◦ A2 + G ;
#        if deg(A) ≤ 1 then
#            Results ← Results ∪ {(A2, A + F(c2))} ;
#        end
#    end
#end
#return Results
