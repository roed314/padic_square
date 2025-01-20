Attach("Montes.m");
Attach("Milstead.m");


Z3 := pAdicRing(3,30);
Z3x<x> := PolynomialRing(Z3);
phi3 := x^9+3*x^8+3*x^6+9*2*x^4+6+9*2;
GaloisGroupMilstead(phi3);

Z3 := pAdicRing(3,30);
U<a> := UnramifiedExtension(Z3,4);
Uy<y> := PolynomialRing(U);
f :=  y^9 + 3*a *y^5 + 6*a;
GaloisGroupMilstead(f);

Z3 := pAdicRing(3,30);
Z3x<x> := PolynomialRing(Z3);
T<a> := TotallyRamifiedExtension(Z3,x^2+3);
Ty<y> := PolynomialRing(T);
f :=  y^9 + 3*a*y^5 + a;
GaloisGroupMilstead(f);

Z2 := pAdicRing(2,300);
Z2x<x> := PolynomialRing(Z2);
f := x^16 + 16 ;
GaloisGroupMilstead(f);

Z2 := pAdicRing(2,30);
Z2x<x> := PolynomialRing(Z2);
f := x^16 + 2*x^15 + 2*x^13 + 2*x^10 + 2*x^8 + 2;
GaloisGroupMilstead(f);
