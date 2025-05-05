// Usage: ls /scratch/lf/to_filter/ | parallel -j100 --timeout 300 --results /scratch/lf/filterout/ "magma -b fname:={0} ComputePolRedPadicNew.m"

Attach("../Pauli/polredabs.m");
SetColumns(0);
SetVerbose("Monge", 2);

function remove_whitespace(X)
    return Join(Split(Join(Split(X," "),""),"\n"),"");
end function;
function sprint(X)
    return remove_whitespace(Sprint(X));
end function;

infile := Sprintf("/scratch/lf/to_filter/" * fname);
outfile := Sprintf("/scratch/lf/filtered/" * fname);
pieces := Split(fname, ".");
p := StringToInteger(pieces[1]);
n := StringToInteger(pieces[2]);
f := StringToInteger(pieces[3]);
ZZ := Integers();
abspoly, relpoly := Explode(Split(Read(infile), "|"));
if f eq n then
    // Do nothing to abspoly, since we already have the conway_or_jr_polynomial
    abspoly := remove_whitespace(abspoly);
    relpoly := "x - " * Sprint(p);
else

    prec := 4*n;

    k0 := pAdicRing(p, prec);
    R0<x> := PolynomialRing(k0);
    abspoly := eval abspoly;
    abspoly := PolRedPadic(abspoly);
    Zx<x> := PolynomialRing(ZZ);
    abspoly := sprint(Zx!abspoly);

    if f eq 1 then
        relpoly := abspoly;
    else
        k<t> := UnramifiedExtension(k0, conway_or_jr_polynomial(k0, f));
        Rrel<x> := PolynomialRing(k);
        relpoly := eval relpoly;
        relpoly := PolRedPadic(relpoly);
        Zt<t> := PolynomialRing(ZZ);
        AssignNames(~Zt, ["t"]);
        Ztx<x> := PolynomialRing(Zt);
        relpoly := sprint(Ztx![Zt![ZZ!a : a in Coefficients(c)] : c in Coefficients(relpoly)]);
    end if;
end if;
PrintFile(outfile, Sprintf("%o|%o", abspoly, relpoly));
quit;
