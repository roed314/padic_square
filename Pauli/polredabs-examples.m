Attach("polredabs.m");

////////////////////////////////////////////////
// PolRedPadic for polynomials over the integers

ZX<X> := PolynomialRing(Integers());         

// Eisenstein over Z2 and Z3
f := X^6+246*X^4+84*X+30; 
PolRedPadic(f,2);        
PolRedPadic(f,3);     

// Eisenstein form over Z3
g  := Evaluate(f,X^2+2*X+2)+3*X;
G := PolRedPadic(g,3); 
G;
// LaTeX output of rg 
String(G:Latex);

////////////////////////////////////////////////
// PolRedPadic for polynomials over local rings

z3 := pAdicRing(3,30);
z3x<x> := PolynomialRing(z3);
f := x^6+6*x^5+18*x+3;
F := PolRedPadic(f);

// relatively reduced
// over an unramified extension
U<a> := UnramifiedExtension(z3,2);
Uy<y> := PolynomialRing(U);
g := y^6+6*y^5+18*y+3*a;
G := PolRedPadic(g);

// invariants
// residual polynomials are as in the LMFDB
ResidualPolynomial(g);  
ResidualPolynomialClasses(g);  
ResidualPolynomialDistinguished(g);
// ramification polygon and polynomial
rp, rho := RamificationPolyAbs(g);
// ramification polygon is reverse of the ramification polygon in the LMFDB
Slopes(rp);
Vertices(rp);

// absolute reduced polynomial
phi,nu,alpha := EisensteinForm(g);
Phi := PolRedPadic(phi);
String(Phi:Latex);

