Attach("Montes.m");
Attach("Milstead.m");
Z3 := pAdicRing(3,100);
Z3x<x> := PolynomialRing(Z3);
phi3 := x^9+3*x^8+3*x^6+9*2*x^4+6+9*2;
GaloisGroupMilstead(phi3);
