# Getting Started

## Installation

- Download a copy of the repository [here](https://github.com/cjdoris/pAdicExtensions).
- In Magma, attach the `spec` file (see the example below).

## Examples

### Main features

```
> // Do this once in the session
> AttachSpec("/path/to/package/spec");
>
> // This package avoids working with explicit fields, and
> // instead uses "template" fields which contain just enough
> // information, such as the prime and ramification degree
> F := TemplatepAdicField(3);
>
> // Find the ramification polygons of degree 9 over Q_3
> // with dicriminant valuation 10+9-1=18
> Ps := AllRamificationPolygons(F, 9 : J0:=10);
> Ps;
[
  Ramification polygon with points [(1, 10)*, (3^2, 0)*],
  Ramification polygon with points [(1, 10)*, (3^1, 3)*, (3^2, 0)*],
  Ramification polygon with points [(1, 10)*, (3^1, 6)*, (3^2, 0)*]
]
>
> // Pick one
> P := Ps[2];
> P;
Ramification polygon with points [(1, 10)*, (3^1, 3)*, (3^2, 0)*]
>
> // Make a template
> T := TemplateForEisensteinPolynomials(P);
> T;
Template for an element of
  Univariate polynomial ring over Template 3-adic field (d=1, f=1, e=1) actually
  3-adic field
with coefficients
   0 (1; F* F)
   1 (2; F*)
   2 (2; F)
   3 (1; F* F)
   4 (2; F)
   5 (2; )
   6 (1; F)
   7 (2; )
   8 (2; )
   9 (0; {1})
> 
> // The number of possible outputs
> #T;
1944
>
> // Generate one of them at random
> // You can also call ToSequence to list all of them
> f := Random(T);
> f;
(0; 1)*x^9 + (1; 1)*x^6 + (2; 2)*x^4 + (1; 1 2)*x^3 + (2; 2)*x + (1; 1 2)
>
> // Coerce it into the actual 3-adic field
> // NOTE: In this case, this is equivalent to Actual(f)
> K := pAdicField(3);
> m := Embedding(Parent(f), PolynomialRing(K));
> m(f);
(1 + O(3^20))*$.1^9 + (3 + O(3^21))*$.1^6 + (2*3^2 + O(3^22))*$.1^4 + (7*3 + O(3^21))*$.1^3 + (2*3^2 + O(3^22))*$.1 + 7*3 + O(3^21)
>
> // Make the corresponding extension
> // NOTE: In this case, this is equivalent to ActualExtension(f)
> L<pi> := Extension(K, f, m);
```

### Completions

```
> // Make a 2-adic completion of the Gaussian rationals
> F<u> := QuadraticField(-1);
> pi := u+1;
> K, e := Completion(F, ideal<Integers(F) | 2, pi>);
>
> // Find out which discriminants are possible
> tK := TemplatepAdicField(K);
> OrePossibilities(tK, 8);
[ 1, 3, 5, 7, 9, 11, 13, 15, ..., 42, 43, 44, 45, 46, 47, 48 ]
>
> // Generate all ramification polygons with discriminant 45, and with all possible residues and CC-residues
> Ps := AllRamificationPolygons(tK, 8 : J0:=45, Res, URes, Classes);
> Ps;
[
  Ramification polygon with points [(1, 45, 1)*, (2^1, 30, 1)*, (2^2, 16, 1)*, (2^3, 0, 1)*] and uniformizer residue 1,
  Ramification polygon with points [(1, 45, 1)*, (2^1, 32, 1)*, (2^2, 16, 1)*, (2^3, 0, 1)*] and uniformizer residue 1
]
>
> // generate templates; check that (in this case) polynomials in the templates correspond 1-1 with extensions
> Ts := [TemplateForEisensteinPolynomials(P) : P in Ps];
> &+[#T : T in Ts], NumberOfExtensions(Integers(K), 8 : F:=1, j:=45);
65536 65536
>
> // generate a random Eisenstein polynomial from one of these
> tf := Random(Ts[1]);
> tf;
(0; 1)*x^8 + (6; 1)*x^7 + (4; 1 1)*x^6 + (6; 1 1)*x^5 + (3; 1 0 0 1 0 1)*x^4 + (8; 1)*x^3 + (1; 1 0 0 0 1)
>
> // Embed the polynomial into the original number field with respect to the uniformizer pi=u+1
> // NOTE: We could have defined tK as TemplatepAdicCompletion(F, ideal<...> : Uniformizer:=pi)
> //   or as TemplatepAdicCompletion(e : Uniformizer:=pi) and then get f := Global(tf)
> te := Embedding(tK, F, e : Uniformizer:=pi);
> f := Polynomial([c@te : c in Coefficients(tf)]);
> f;
x^8 - 8*u*x^7 + (-4*u - 8)*x^6 + (-16*u + 8)*x^5 + (-6*u + 14)*x^4 + 16*x^3 - 3*u - 3
```