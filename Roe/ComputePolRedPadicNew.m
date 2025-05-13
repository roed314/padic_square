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
relpoly := Read(infile);
if "|" in relpoly then
    abspoly, relpoly := Explode(Split(Read(infile), "|"));
elif f eq 1 then
    abspoly := relpoly;
else
    abspoly := ""; // set below
end if;
if f eq n then
    // Do nothing to abspoly, since we already have the conway_or_jr_polynomial
    abspoly := remove_whitespace(abspoly);
    relpoly := "x-" * Sprint(p);
else

    prec := 4*n;
    k0 := pAdicRing(p, prec);
    Zt := PolynomialRing(ZZ);

    if f ne 1 then
        k<t> := UnramifiedExtension(k0, conway_or_jr_polynomial(k0, f));
        Rrel<x> := PolynomialRing(k);
        relpoly := eval relpoly;
    end if;
    if abspoly cmpeq "" then // implies f != 1
        B<b> := TotallyRamifiedExtension(k, relpoly);
        Bx<x> := PolynomialRing(B);
        y := Roots(Bx!DefiningPolynomial(k) - b, B)[1][1];
        abspoly := MinimalPolynomial(y, k0);
    else
        R0<x> := PolynomialRing(k0);
        abspoly := eval abspoly;
    end if;
    abspoly := PolRedPadic(abspoly);
    AssignNames(~Zt, ["x"]);
    abspoly := sprint(Zt!abspoly);
    if f eq 1 then
        relpoly := abspoly;
    else
        relpoly := PolRedPadic(relpoly);
        AssignNames(~Zt, ["t"]);
        Ztx<x> := PolynomialRing(Zt);
        relpoly := sprint(Ztx![Zt![ZZ!a : a in Coefficients(c)] : c in Coefficients(relpoly)]);
    end if;
end if;
PrintFile(outfile, Sprintf("%o|%o", remove_whitespace(abspoly), remove_whitespace(relpoly)));
quit;
