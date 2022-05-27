
/*
unramified a polynomial in t
eisen is a poly in y with coeff containing t
*/

intrinsic myii(ent) -> .
 {Doc}
    p := ent[2];
    uncoef := ent[5];
    eis := ent[4];
    zp:=pAdicRing(p,300);
    uncoef2:= [zp! z : z in uncoef];
    zptt<tt>:=PolynomialRing(zp);
    unpol:=elt<zptt | uncoef2>;
    unr<t>:=UnramifiedExtension(zp, unpol);
    unry<yy>:=PolynomialRing(unr);
    eiscoef := Coefficients(eis);
    eiscoef2 := <Coefficients(z) : z in eiscoef>;
    eiscoef2 := <<zp ! z : z in u> : u in eiscoef2>;
    eiscoef2;
    eiscoef3 := < elt<zptt | z> : z in eiscoef2>;
    eiscoef3;
    eiscoef := [unr! z : z in eiscoef3];
    eis2 := elt<unry | eiscoef>;

    return 1;
end intrinsic;


