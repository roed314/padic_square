// Usage: ls /scratch/lf/filtered | parallel -j80 -- timeout 300 --results /scratch/lf/adout/ "magma -b fname:={0} AutDisc.m"

Attach("../Pauli/polredabs.m");
SetColumns(0);

infile := "/scratch/lf/filtered/" * fname;
outfile := "/scratch/lf/outdisc/" * fname;

s := Read(infile);
abspoly, relpoly := Explode(Split(s, "|"));

p, n, f, i := Explode([StringToInteger(c) : c in Split(fname, ".")]);
prec := 4*n;
K0 := pAdicRing(p, prec);
K := UnramifiedExtension(K0, conway_or_jr_polynomial(K0, f));
Rrel<x> := PolynomialRing(K);
relpoly := eval relpoly;
L := TotallyRamifiedExtension(K, relpoly);
c := Valuation(Discriminant(L, K0));
A := AutomorphismGroup(L, K0);
PrintFile(outfile, Sprintf("%o|%o", c, #A));
quit;
