https://github.com/yue-zhou-ovgu/APN-functions/blob/main/MAGMA_Proofs.txt


This file provides two MAGMA codes for this paper 
[Kangquan Li, Yue Zhou, Chunlei Li, and Longjiang Qu, Two new families of quadratic APN functions, 
accepted by IEEE Transactions on Information Theory, 2022]
		
A. The MAGMA code in the proof of Theorem 2:
		
P<x,y,a,b> := PolynomialRing(GF(2),4);
F := a*x^2+(a^2+b^2+b)*x+(a+b)*y^2+(a+b^2)*y;
G := (a+b)*x^4+b^2*x^2+(a^4+b)*x+b*y^4+a^2*y^2+(a^4+a+b^4)*y;
H := Resultant(F,G,y);
Factorization(H);
				
	
		
		
		

		
B. The MAGMA code in the proof of Lemma 5.

In the equation system (11),  we assume that $a1 := a^{2^m}$, 
$a2 := a^{2^{2m}}$, $as := a^{2^s}$, $as1 := a^{2^{m+s}}$, $as2 := a^{2^{2m+s}}$,
$u := \mu$, $u1 := \mu^{2^m}$ and $u2 := \mu^{2^{2m}}$. 
Then $A, B, C, D, E$ can be seen as polynomials in $\gf_2[a,a1,a2,as,as1,as2,u,u1,u2,v]$. 
We also suppose $X1 := X^{2^m}$ and  $X2 := X^{2^{2m}}$ for $X\in\{A,B,C,D,E\}$. 
Then we can use the following MAGMA code to compute the factorizations of $U1,U2,U3,V1,V2,V3$. 
		
	
P<a,a1,a2,as,as1,as2,u,u1,u2,v> := PolynomialRing(GF(2),10); 
A := (as1+u*as+a)*as2;
A1 := (as2+u1*as1+a1)*as;
A2 := (as+u2*as2+a2)*as1;
B := (as2+u1*as1+a1+u1*(as1+u*as+a))*as1;
B1 := (as+u2*as2+a2+u2*(as2+u1*as1+a1))*as2;
B2 := (as1+u*as+a+u*(as+u2*as2+a2))*as;
C := (as1+u*as+a+v*a)*a1;
C1 := (as2+u1*as1+a1+v*a1)*a2;
C2 := (as+u2*as2+a2+v*a2)*a;
D := u*(as2+u1*as1+a1)*as;
D1 := u1*(as+u2*as2+a2)*as1;
D2 := u2*(as1+u*as+a)*as2;
E := (as2+u1*as1+a1+v*a1)*a;
E1 := (as+u2*as2+a2+v*a2)*a1;
E2 := (as1+u*as+a+v*a)*a2;
			
U1 := D2*E1*E+A*C2*E1+B1*C2*C;
U2 := A2*E1*E+B*C2*E1+C2*C*D1;
U3 := B2*E1*E+C2*D*E1+A1*C2*C;
V1 := A2*A^2*C1+A*B*C1*D2+A*B1*B*E2+A^2*D1*E2;
V2 := A2*A^2*E1+A*B*D2*E1+A2*A*B1*C+A*C*D2*D1;
V3 := A2*A*B1*E+A*B1*B*C2+A^2*C2*D1+A*D2*D1*E;
			
Factorization(U1);
Factorization(U2);
Factorization(U3);
Factorization(V1);
Factorization(V2);
Factorization(V3);