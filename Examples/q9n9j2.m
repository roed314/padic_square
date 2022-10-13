//Loading "load.m"
z3 := pAdicRing(3,20);
z9<a> := UnramifiedExtension(z3,2);
z9x<x> := PolynomialRing(z9);
// K = Z3(a) = Z3[x] / x^2 + 2*x + 2
field := recformat <K,n,j,I,R,A,f,aut,gal>;
L := [
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + 1>, f:=x^9 + 3*x^2 + 3, aut:="2T1", gal:="9T9">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a>, f:=x^9 + 3*a*x^2 + 3, aut:="1T1", gal:="9T15">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a^2>, f:=x^9 + (3*a + 3)*x^2 + 3, aut:="1T1", gal:="9T9">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a^3>, f:=x^9 + (6*a + 3)*x^2 + 3, aut:="1T1", gal:="9T15">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + 2>, f:=x^9 + 6*x^2 + 3, aut:="2T1", gal:="9T9">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a^5>, f:=x^9 + 6*a*x^2 + 3, aut:="1T1", gal:="9T15">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a^6>, f:=x^9 + (6*a + 6)*x^2 + 3, aut:="1T1", gal:="9T9">,
rec< field | K:="Q9", n:=9, j:=2, I:=[ 2, 2, 0 ], R:=[ <1, 2>, <9, 0> ], A:=<x^2 + a^7>, f:=x^9 + (3*a + 6)*x^2 + 3, aut:="1T1", gal:="9T15">,
rec< field |>
];
