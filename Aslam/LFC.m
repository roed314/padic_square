


//import "Module.m" : lattice;
//import "NormEquation" : ClNormEquation;
import "brauer.m"  : compute_LPmul_modX;
//////////////////////////CLocalFundamentalClassSerre//////////////


// verbose comments
declare verbose Brauer, 2;
declare verbose CocycleLFC, 2;
//SetVerbose("CocycleLFC", 0);
declare verbose FrobEq, 3;
//SetVerbose("FrobEq", 1);

/*
intrinsic automorphismGroup(OL::RngPad, OK::RngPad) -> GrpPerm, Map, SeqEnum
{ Compute the automorphism group of a local field extension;
  similar to the Magma function AutomorphismGroup, but also
  works, if the fields are of infinite precision. }
    return automorphism_group(OL, OK);
end intrinsic;


intrinsic automorphismGroup(OL::FldPad, OK::FldPad) -> GrpPerm, Map, SeqEnum
{ Compute the automorphism group of a local field extension;
  similar to the Magma function AutomorphismGroup, but also
  works, if the fields are of infinite precision. }
    return automorphism_group(OL, OK);
end intrinsic;
// See automorphismGroup.
 
function automorphism_group(OL, OK)
    local p, f, rts, r, Aut, permut, i, g, G, H, x, psi;

    p, f := primitive_element(OL, OK);

    rts := Roots(f, OL);
    rts := [r[1] : r in rts];
    //assert #rts eq Degree(OL, OK);

    Aut := Automorphisms(OL, OK);

    //G := PermutationGroup< {1..Degree(OL,OK)} | [[Index(rts,r) : r in  [g(r) : r in rts] ] : g in Aut] >;
    permut := [   [   [i : i in [1..#rts] | Valuation(rts[i]-r) eq Minimum(Precision(r),Precision(rts[i])) ][1]
             : r in  [g(r) : r in rts] ]
         : g in Aut ];
    G := PermutationGroup< {1..#rts} | permut >;
    H := [G!x : x in permut];

    psi := map< G -> Aut | x :-> Aut[Index(H,x)] >;

    return G, psi, Aut;
end function;

*/



function ali_cocycle_lfc_G(L, K, psi, precision, NoNorm)
    local steps, G, OL, Zp, Qp, pi, pi_sigma, pisigmapi, g, i, u_sigma, phi,
          e, OL2ub, OL2, OL1, L1, K1, incl, GG, AutL, AutL1, sigma, psi1,
          OLmal, m, u, bool, gamma1;
if psi cmpeq 0 then
        _,psi,_ := AutomorphismGroup(L,K);
        psi := map<Domain(psi) -> Codomain(psi) | g :-> psi(g^(-1))>;
    end if;
is_equal := func< x,y | ChangePrecision(x,m) eq ChangePrecision(y,m)
                        where m := Minimum(Precision(x), Precision(y)) >;


steps := precision+2;
    G := Domain(psi);
    OL := RingOfIntegers(L);
    Zp := pAdicRing(RingOfIntegers(OL));
    Qp := FieldOfFractions(Zp);
    d := InertiaDegree(OL,Zp);


if d eq Degree(L,K) then
        pi := UniformizingElement(K);
        phi := FrobeniusAutomorphism(L,K);
	//phi := FrobeniusAutomorphism(L);
        phi := [g : g in G | &and([ is_equal(psi(g)(b), phi(b)) : b in generators_pad(L,K) ]) ];
        if #phi eq 0 then
            // ungenauer
            phi := FrobeniusAutomorphism(L, K);
            phi := [g : g in G | &and([ Valuation(psi(g)(b) - phi(b)) ge L`DefaultPrecision*95/100 : b in generators_pad(L,K) ]) ];
            //print d;
            //print FrobeniusAutomorphism(L, K)(L.1);
            //print [psi(g)(L.1) : g in G ];
            //print G;
        end if;
        phi := phi[1];
        GG := [phi^i : i in [0..#G-1]];
        return map< car<G,G> -> L | x :-> (Index(GG, x[1])+Index(GG,x[2])-2 lt #G select L!1 else pi) >;
    end if;

// is_equal and generators_pad;


if steps gt L`DefaultPrecision then
        error "Precision of p-adic field L not high enough to compute the cocycle!";
    end if;
    if steps gt Minimum([ Precision(r[1]) : r in Roots(DefiningPolynomial(L),L)]) then
        error "Precision of automorphisms of L not high enough to compute the cocyle!";
    end if;

 if NoNorm then
        pi := UniformizingElement(L);
        // Frobenius-Gleichungen loesen
        pi_sigma  := [ OL!psi(g)(pi) : g in G];
        pisigmapi := [ OL!(pi_sigma[i]/pi) : i in [1..#pi_sigma]];
        vprint CocycleLFC, 1: "Solve Frobenius equations without norms";
        vtime  CocycleLFC, 1: u_sigma, phi := FrobeniusEquation(pisigmapi, steps);

        // Erweiterungskoerper
        OL2 := Parent(u_sigma[1]);
        // Galois-Action on \prod OL2
        frobAction, GAction, frobIndex := galois_act_L_tensor_Knr(OL, OL2, psi, phi);
else

        if Precision(L) eq Infinity() then
            e := RamificationIndex(L,K);
            OL2ub := UnramifiedExtension(OL,e);
            vprint CocycleLFC, 1: "Switch to bounded precision, prec = ", OL2ub`DefaultPrecision;
            OL2 := ChangePrecision(OL2ub, OL2ub`DefaultPrecision);
            OL1 := BaseRing(OL2);
            L1 := FieldOfFractions(OL1);
            //L1 := ChangePrecision(L, L`DefaultPrecision);
            //K1 := BaseField(L1);
            //K1 := ChangePrecision(K, K`DefaultPrecision);
            K1 := pAdicField(L1);

            incl := Coercion(L,L1);
            G := Domain(psi);
            GG := [g : g in G];
            AutL := [psi(g) : g in GG];
            AutL1 := [ (incl^(-1))*sigma*incl  : sigma in AutL];
            psi1 := map< G -> AutL1 | g :-> AutL1[Index(GG,g)] >;

            if steps gt Minimum([ Precision(r[1]) : r in Roots(DefiningPolynomial(L1),L1)]) then
            //if steps gt Minimum([ Precision(psi1(g)(L1.1)) : g in G]) then
                error "Precision of automorphisms of L (bounded) not high enough to compute the cocyle!";
            end if;
 else
            L1 := L;
            K1 := K;
            e := RamificationIndex(L1,K1);
            OL1 := RingOfIntegers(L1);
            OL2 := ext<OL1 |e>;
        end if;
    
     vprintf CocycleLFC, 1: "Compute UnitGroup... ";
       // vtime   CocycleLFC, 1: OLmal, m := UnitGroup(OL1);
     //        Attach("NormEquation.m");
        u := UniformizingElement(K1)/UniformizingElement(L1)^e;
        vprintf CocycleLFC, 1: "Solve Norm equation... ";
       // vtime   CocycleLFC, 1: bool, gamma1 := MyNormEquation(OL2,m,OL1!u);
        vtime   CocycleLFC, 1: gamma1 := ClNormEquation(OL2,OL1!u);   
       
  /*      vprintf CocycleLFC, 1: "Compute UnitGroup... ";
        vtime   CocycleLFC, 1: OLmal, m := UnitGroup(OL1);

        u := UniformizingElement(K1)/UniformizingElement(L1)^e;
        vprintf CocycleLFC, 1: "Solve Norm equation... ";
        vtime   CocycleLFC, 1: bool, gamma1 := NormEquation(OL2,m,OL1!u);
*/
        //assert bool;
        pi:= gamma1*UniformizingElement(L1);
if Precision(L) eq Infinity() then
            vprint CocycleLFC, 1: "Switch back to unbounded precision";
            pi := OL2ub!pi;
            OL2 := OL2ub;
        else
            pi := OL2!pi;
        end if;
       phi := FrobeniusAutomorphism(OL2, OL);
 //    phi := FrobeniusAutomorphism(OL2);// "this makes problem in the following action";
        frobAction, GAction, frobIndex := galois_act_L_tensor_Knr(OL, OL2, psi, phi);
        pi_sigma := [GAction(g,<pi : i in [1..d]>)[1] : g in G];
        pisigmapi := [ OL2!(pi_sigma[i]/pi) : i in [1..#pi_sigma]];
        vprintf CocycleLFC, 1: "Solve Frobenius equations... ";
        vtime  CocycleLFC, 1: u_sigma, phi := FrobeniusEquation(pisigmapi, steps, OL);
    end if;

 if GetVerbose("CocycleLFC") ge 2 then
        vprint CocycleLFC, 2: "Test FrobeniusEquation result";
        assert &and({Valuation(phi(u_sigma[i])/u_sigma[i] - pisigmapi[i]) ge steps : i in [1..#u_sigma]});
    end if;

    // Kozykel
    d := InertiaDegree(OL,Zp);
    L2 := FieldOfFractions(OL2);
    prodL2 := Domain(frobAction);
    prodOL2 := CartesianProduct([OL2 : y in [1..d] ]);

    tup_sigma := [];
    GG := [g : g in G];
    for g in GG do
        ind := Index(GG,g);
        frobIdx := frobIndex[ind];
        if frobIdx eq 0 then
            frobIdx := d;
        end if;
        Append(~tup_sigma, prodOL2! < u_sigma[ind] * ( i le frobIdx select 1 else  pi_sigma[ind] ) : i in [1..d] > );
    end for;

 c := map<car<G,G> -> prodL2   | x  :->   // x = <x[1], x[2]>
                tupelQuotient(
                    tupelProduct(
                        GAction(x[1], tup_sigma[Index(GG, x[2])]),
                        tup_sigma[Index(GG, x[1])]
                    ),
                    tup_sigma[Index(GG, x[1]*x[2])]
                )
            >;

 c := precompute_map(c);
    // Komponenten gleich modulo pi^(precision+1)
    assert Minimum([ Minimum([Valuation(y[1]-y[i]) : i in [1..#y]])
        where y is c(x) :  x in Domain(c)]) ge (precision+1);
    if Degree(Codomain(c)[1]) ge 2 then
        // erste Komponente in L modulo pi^(precision+1)
        assert Minimum([ Minimum([ Valuation(z) : z in ElementToSequence(y[1])[2..Degree(Parent(y[1]))]])
            where y is c(x)  :  x in Domain(c)]) ge (precision+1);
    end if;
 gamma := map< Domain(c) -> FieldOfFractions(L) | x :->  ( elem_to_seq(c(x)[1], L)[1] )^(-1) >;

    return gamma;
end function;


function restrictAutomorphism(sigma, K)
    incl := Coercion(K, Domain(sigma));
    return hom<K-> K | x:-> sigma(x), y:-> (sigma^(-1))(y)>;
end function;


is_zero  := func< x   | x eq ChangePrecision(Zero(Parent(x)),Precision(x)) >;
is_equal := func< x,y | ChangePrecision(x,m) eq ChangePrecision(y,m)
                        where m := Minimum(Precision(x), Precision(y)) >;


intrinsic generators_pad(L, K)->.{}
    G := [*  *];
    E := L;
    while not E eq K do
        if E.1 eq Zero(E) then
            Append(~G, One(E));
        else
            Append(~G, E.1);
        end if;
        E := BaseRing(E);
        //Insert(~G,1, E.1);
    end while;
    return G;
end intrinsic;

/*
intrinsic automorphismGroup(OL::RngPad, OK::RngPad) -> GrpPerm, Map, SeqEnum
{ Compute the automorphism group of a local field extension;
  similar to the Magma function AutomorphismGroup, but also
  works, if the fields are of infinite precision. }
    return automorphism_group(OL, OK);
end intrinsic;
*/

/* See automorphismGroup.
 */
/*
intrinsic automorphism_group(OL, OK)->.{}
    local p, f, rts, r, Aut, permut, i, g, G, H, x, psi;

    p, f := primitive_element(OL, OK);

    rts := Roots(f, OL);
    rts := [r[1] : r in rts];
    //assert #rts eq Degree(OL, OK);

    Aut := Automorphisms(OL, OK);

    //G := PermutationGroup< {1..Degree(OL,OK)} | [[Index(rts,r) : r in  [g(r) : r in rts] ] : g in Aut] >;
    permut := [   [   [i : i in [1..#rts] | Valuation(rts[i]-r) eq Minimum(Precision(r),Precision(rts[i])) ][1]
             : r in  [g(r) : r in rts] ]
         : g in Aut ];
    G := PermutationGroup< {1..#rts} | permut >;
    H := [G!x : x in permut];

    psi := map< G -> Aut | x :-> Aut[Index(H,x)] >;

    return G, psi, Aut;
end intrinsic;

*/
///changed myself
/*
function frobenius_automorphism(L, K)
    phi := map<K->K | x :-> x>;
    if L eq K then
        return phi;
    else
        q := #ResidueClassField(RingOfIntegers(K));
        return continue_frobenius(phi, L, q);
    end if;
end function;
*/

/* 
 * Given the automorphism map psi: G-> Aut(L) of a p-adic extension
 * L/K. Compute a group element that represents the Frobenius
 * automorphism. 
 */

/////////////////////////////changed////////////////

function frobenius_automorphism_grpelt(psi)
    local G, L, K, B;
    G := Domain(psi);
    L := Domain(psi(G!1));
  //B:=[**];
    B := generators_pad(L,K);
    // find base field K of the automorphism group and 
    // generate list of generating elements of L/K.
    K := BaseRing(L);
  Qp := PrimeRing(L);
  /* if Type(L) eq RngPad then
        Qp := pAdicRing(L);
    else
        Qp := pAdicField(L);
    end if;
   */
/*
while Degree(L,K) lt #G do
        Append(~B, K.1);
        K := BaseRing(K);
        if K cmpeq Qp and Degree(L,K) lt #G then
            error "Couldn't find fixed field of automorphism group in field tower of L/Qp!";
        end if;
    end while;
*/
    OK := RingOfIntegers(K);
    q := #ResidueClassField(OK);

    OL := RingOfIntegers(L);
    R := ResidueClassField(OL);

    for g in G do
        sigma := psi(g);
        if &and([ R!sigma(b) - (R!b)^q eq R!0 : b in B  ]) then
            return g;
        end if;
    end for;

    error "No Frobenius element found!";
end function;
/*
intrinsic primitive_element(L, K)->.{}
	//print "primitive_element";
    E := BaseRing(L);
    if E eq K then
        return L.1, MinimalPolynomial(L.1);
    else
        alpha := L.1;
        //print("Berechne primitives Element rekursiv");
        beta := primitive_element(E,K);
        //print beta;
        n := Degree(L,K);
        for j in [1..10] do
            gamma := Random(1,10)*alpha + beta;
            //print "berechne MiPo";
            f := MinimalPolynomial(gamma,K);
            if Type(K) eq RngPad then
                f := PolynomialRing(FieldOfFractions(K))!f;
            end if;
            try
                //print "faktorisiere quadratfrei";
                //f := SquareFreeFactorization(f)[1,1];
                //print "faktorisiere";
                if Degree(Factorization(f)[1,1]) eq n then
                    return gamma, f;
                end if;
            catch e
                gamma := 0;
            end try;
        end for;
        error "Did not find a primitive element!";
    end if;
end intrinsic;

*/

//frobenius_automorphism related

intrinsic continue_frobenius(frob, L, q)->.{}

    local K, B, b, f, rts, cont, x, r, m, i, cand;
    //require (Type(L) eq FldPad) or 
    //        (Type(L) eq RngPad): 
    //        "Bad argument types\np-Adic Rings or Fields required";
    
    
    K := Domain(frob);
    
    B := Reverse(generators_pad(L,K));
    for b in B do
        // continuations of frob with values in L
        f := MinimalPolynomial(b);
	// "ali_computation of roots is expensive  try with Magma'FrobAut";
        rts := {r[1] : r in Roots( Polynomial([ frob(x) :  x in Coefficients(f)]), L) };
//"Ali Changed";	
        assert #rts eq Degree(f);
        cont := [ map< Parent(b) -> L | 
            x :-> &+([ frob(elem[i])*r^(i-1) : i in [1..#elem]] 
            where elem is elem_to_seq(x,BaseRing(Parent(b)))) >  
            : r in rts 
        ];
        
//where elem is elem_to_seq(x,BaseRing(Parent(b)))) >	
        for m in cont do
            cand := true;
            for i in [0..Degree(Parent(b))-1] do
                if Valuation(m(b^i)-(b^i)^q) eq 0 then
                    cand := false;
                    break;
                end if;
            end for;
            if cand then
                // found a lift of the frobenius
                frob := m;
                break;
            end if;
        end for;
        
        if not cand then
            error "No Frobenius found!";
        end if;
    end for;
    
    return frob;
end intrinsic;


/*
intrinsic FrobeniusAutomorphism(L::FldPad, K::FldPad) -> Map
{ Computes the Frobenius automorphism of K, restricted to L/K.}
    return frobenius_automorphism(L,K);
end intrinsic;

intrinsic FrobeniusAutomorphism(L::RngPad, K::RngPad) -> Map
{ Computes the Frobenius automorphism of K, restricted to L/K.}
    return frobenius_automorphism(L,K);
end intrinsic;
*/
intrinsic FrobeniusAutomorphism(psi::Map) -> GrpElt
{ Given the automorphism map psi: G-> Aut(L) of an extension
  L/K. Compute the group element that represents the Frobenius
  automorphism. }
    return frobenius_automorphism_grpelt(psi);
end intrinsic;


intrinsic FrobeniusEquation(c::RngPadElt, precision::RngIntElt) -> RngPadElt, Map
{Solves the equation x^(phi-1)=c, c in OE^\times, up to the given precision,
where phi is the Frobenius automorphism of OK.
The solution x and the automorphism phi are returned.
If a sequence C of elements is given, a sequence of solutions is returned.
If OK is not given, OK=OE. Otherwise, OE must be an extension of OK.
Note, that whenever the norm of c over OK is not 1, this can generate
huge extensions of OE.}

    local x, frob;
    x, frob := frobenius_equation(c, precision);
    return x, frob;
end intrinsic;

intrinsic FrobeniusEquation(C::SeqEnum[RngPadElt], precision::RngIntElt) -> SeqEnum[RngPadElt], Map
{Solves the equation x^(phi-1)=c, c in OE^\times, up to the given precision,
where phi is the Frobenius automorphism of OK.
The solution x and the automorphism phi are returned.
If a sequence C of elements is given, a sequence of solutions is returned.
If OK is not given, OK=OE. Otherwise, OE must be an extension of OK.
Note, that whenever the norm of c over OK is not 1, this can generate
huge extensions of OE.}

    local x, frob, c, X;

    c := C[1];
    q := #ResidueClassField(Parent(c));
    x, frob := frobenius_equation(c, precision);
    X := [x];

    for i in [2..#C] do
        O := Parent(x);
        c := C[i];
        x, frob := frobenius_equation_OK(O!c, precision, frob, q);
        X := X cat [x];
    end for;

    return X, frob;
end intrinsic;


intrinsic FrobeniusEquation(c::RngPadElt, precision::RngIntElt, OK::RngPad) -> RngPadElt, Map
{Solves the equation x^(phi-1)=c, c in OE^\times, up to the given precision,
where phi is the Frobenius automorphism of OK.
The solution x and the automorphism phi are returned.
If a sequence C of elements is given, a sequence of solutions is returned.
If OK is not given, OK=OE. Otherwise, OE must be an extension of OK.
Note, that whenever the norm of c over OK is not 1, this can generate
huge extensions of OE.}

    local q, frob, x;

    q := #ResidueClassField(OK);
    frob :=  map<OK -> OK | x :-> x>;
    frob := continue_frobenius(frob, Parent(c), q);
    x, frob := frobenius_equation_OK(c, precision, frob, q);
    return x, frob;
end intrinsic;

intrinsic FrobeniusEquation(C::SeqEnum[RngPadElt], precision::RngIntElt, OK::RngPad) -> SeqEnum[RngPadElt], Map
{Solves the equation x^(phi-1)=c, c in OE^\times, up to the given precision,
where phi is the Frobenius automorphism of OK.
The solution x and the automorphism phi are returned.
If a sequence C of elements is given, a sequence of solutions is returned.
If OK is not given, OK=OE. Otherwise, OE must be an extension of OK.
Note, that whenever the norm of c over OK is not 1, this can generate
huge extensions of OE.}

    local q, x, frob, c, X;

    q := #ResidueClassField(OK);
    frob :=  map<OK -> OK | x :-> x>;
    frob := continue_frobenius(frob, Parent(C[1]), q);

    c := C[1];
    x, frob := frobenius_equation_OK(c, precision, frob, q);
    X := [x];

    for i in [2..#C] do
        x, frob := frobenius_equation_OK(Parent(x)!C[i], precision, frob, q);
        X := X cat [x];
    end for;

    return X, frob;
end intrinsic;

intrinsic frobenius_equation(c, precision)->.{}
    OK := Parent(c);
    frob := map<OK -> OK | x :-> x>;
    q := #ResidueClassField(OK);
    x, frob := frobenius_equation_OK(c, precision, frob, q);
    return x, frob;
end intrinsic;


/* 
 * Given c in OE, E/K unramified, the Frobenius phi of \tilde K/K and 
 * the characteristic q of the residue class field of K.
 * Solve x^{phi-1}=c up to the given precision.
 */
intrinsic frobenius_equation_OK(c, precision, frob, q)->.{}

    O := Domain(frob);
    pi := UniformizingElement(O);
    F := ResidueClassField(O);
    R<x> := PolynomialRing(F);
    pp := Factorization(x^q-F!c*x);
    f := pp[2,1];

    if Degree(f) gt 1 then
        vprintf FrobEq, 1: "-%o-", Degree(f);
        // Erzeuge Erweiterung in der f eine Nullstelle hat
        O  := ext<O | Degree(f)>;
        pi := UniformizingElement(O);
        F  := ResidueClassField(O);
    else
        //print "Keine Koerpererweiterung notwendig";
    end if;

    R<x> := PolynomialRing(O);
    rts := Roots(R!f);
    x := rts[1,1];
    if Degree(f) gt 1 then
        frob := continue_frobenius(frob, O, q);
    end if;

    // Weitere Schritte
    x, frob := frobenius_equation_cont(c, precision, x, 1, frob, q);

    return x, frob;
end intrinsic;


intrinsic frobenius_equation_cont(c, precision, x, r, frob, q)->.{}

    xi := [x];
    O := Parent(x);
    pi := UniformizingElement(O);
    F := ResidueClassField(O);

    // initialisiere Werte a,b fuer den naechsten Schritt
    a := O!(c*x/frob(x));
    bool, b := IsExactlyDivisible(a-1, pi); // ^ r  !!!
    if not bool then
        error "nicht teilbar 1!";
    end if;

    // weitere Schritte
    for j in [r+1..precision] do
        R<X> := PolynomialRing(F);
     //"Factorisation below consumes enogh time and memory so work on vector space";
	pp := Factorization(X^q-X-F!O!b);
        f := pp[1,1];

        // Falls noetig, erzeuge Erweiterung
        if Degree(f) eq 1 then
            //print "Keine Erweiterung notwendig in Schritt", i;
        else
            //print "Erweiterung vom Grad", Degree(f);
            vprintf FrobEq, 1: "-%o-", Degree(f);
            O := ext<O | Degree(f)>;
        end if;
        pi := UniformizingElement(O);
        F  := ResidueClassField(O);
        rts := Roots(f, O);
        y := rts[1,1];

        lastfrob := frob;
        if Degree(f) gt 1 then
            frob := continue_frobenius(frob, O, q);
        end if;

        z := 1+y*pi^(j-1);
	 // if not (Valuation(frob(z)/z-a[i]) ge #xi[i]+1) then
        //     print "\nFehler: Die Loesung in Schritt", i, "erfuellt die Valuation-Bedingung nicht!\nTeilergebnis wird ausgegeben.\n";
        //     return [ &*(xi[i])  : i in [1..#C]], lastfrob;
        // end if;
        xi := (xi cat [z]);

        if j lt precision then
            // initialisiere Werte a,b fuer naechsten Schritt
            a := (O!a) * z/frob(z);
            bool,bb:=IsExactlyDivisible(O!(a-1),pi^j);
            if bool then
                b := O!bb;
            else
                error "nicht teilbar 2!";
            end if;
        end if;
    end for;

    return &*(xi), frob;
end intrinsic;


intrinsic galois_act_L_tensor_Knr(OL, OL2, psi, phi)->.{}
    local OK, Zp, G, GG, g, d, sigmaHut, frobIndex, L2, prodL2,
        frobeniusMap, Gaction, x;


    // Test whether algorithm is applicable
    OK := maximal_unram_subext_simple(OL);
  //  OK := BaseRing(OL);
    Zp := BaseRing(OK);

    if not generates_maximal_unram(OL2, OL, OK) then
        error "The maximal unramified extension in OL2/Zp cannot be deduced from OL2/OL!";
    end if;

    // Initialization
    G := Domain(psi);
    GG := [g : g in G];
    d := InertiaDegree(OL,Zp);

    // Compute i such that sigma_OK = phi^i for all sigma in G
    // and extensions sigmaHut of sigma such that sigmaHut^(-1)=phi^(-i) on K2
    sigmaHut, frobIndex := continuations_with_unram_restriction(OL, G, psi, OL2);

    // Frobenius automorphism on \prod L2
    L2 := FieldOfFractions(OL2);
    prodL2 := CartesianProduct([L2 : y in [1..d] ]);
    frobeniusMap := map< prodL2 -> prodL2 | x :->  < i eq d select phi(x[1]) else x[i+1]  : i in [1..d]> >;

    // action of G on \prod L2
    Gaction := map< car<G, prodL2> -> prodL2 | x :->
        apply_map( (frobeniusMap^((d-frobIndex[Index(GG,x[1])]) mod d))( x[2] ), sigmaHut[Index(GG, x[1])]) >;

    return frobeniusMap, Gaction, frobIndex;
end intrinsic;

intrinsic apply_map(y, m)->.{}
    local n, i;
    
    case Type(y):
        when Tup:
            n := NumberOfComponents(Parent(y));
        when SeqEnum:
            n := #y;
    else:
      error "Typ nicht unterstuetzt.";
    end case;
    
    for i in [1..n] do
        y[i]:= m(y[i]);
    end for;
    return y;
end intrinsic;



intrinsic continuations_with_unram_restriction(OL, G, psi_OL_Zp, OL2)->.{}
    local Hom_OL_Zp, OK, Zp, sigma, sig, f,
          sigmaHut, frobIndex, sigmaKnrExponentInv;
    
    //require Type(Domain(psi_OL_Zp(G.1))) in {RngPad, FldPad} :
    //        "Bad argument types\nAutomorphisms of p-adic rings/fields needed.";
    
    Hom_OL_Zp := Codomain(psi_OL_Zp);
    OK := maximal_unram_subext_simple(OL);
   //OK := BaseRing(OL);
    Zp := BaseRing(OK);
    
    if not generates_maximal_unram(OL2, OL, OK) then
        print "ACHTUNG: Fuer den angegebenen Koerperturm ist der schnelle Algorithmus nicht anwendbar!";
        print "Benutze den langsamen Algorithmus via primitiven Elementen.";
        
        error "Could not compute 'unramified behaviour' of extension";
    end if;
    
    // Jetzt ist folgendes bekannt:
    // Seien E_i die Zwischenkoerper von OL2/OL mit E_0=OL und E_n=OL2.
    // Dann wird die maximal unverzweigte Teilerweiterung Knr von OL2/OK
    // erzeugt von {E_1.1, E_2.1,...,E_n.1} und Knr(E_0.1)=OL2.
    
    // die Automorphismen von OL/Zp
    sigma := [psi_OL_Zp(g) : g in G];
    if Type(Domain(psi_OL_Zp(G.1))) eq FldPad then
        // Automorphismen auf Ganzheitsring einschraenken
        sigma := [map<OL -> OL | x :-> sig(x) > : sig in sigma];
    end if;
    
    d := Degree(OK);
    d2 := InertiaDegree(OL2, Zp);
    
    if d eq 1 then
        // Ausgangssituation voll verzweigt
        
        frobIndex := [Zero(Integers()) : sig in sigma];
        sigmaKnrExponentInv := frobIndex;
        
        // setze sigma fort mit Identitaet auf unverzweigten Teil
        sigmaHut := [];
        for i in [1..#sigma] do
            sig := sigma[i];
            B := generators_pad(OL2,OL);
            for r in Reverse(B) do
                sig := map<Parent(r) -> Parent(r) | x:-> continued_automorphism_image(x, sig, r) >;
            end for;
            Append(~sigmaHut, sig);
        end for;
        
    elif OL eq OK then
        // Ausgangssitugation unverzweigt
        f := FrobeniusAutomorphism(OK);

       // f := FrobeniusAutomorphism(OK, Zp);
        frobIndex := [find_power(sig, f, OK.1 , InertiaDegree(OK, pAdicRing(OK))) : sig in sigma];
        fInv := inverseAutomorphism(f);
        
        sigmaHut := sigma;
        sigmaKnrExponentInv := [(d-i) mod d : i in frobIndex];
        
    else
        // Verhalten von sigma\in\Aut(OL,Zp) auf OK bestimmen
       // f := FrobeniusAutomorphism(OK, Zp);
       f := FrobeniusAutomorphism(OK);
       frobIndex := [find_power(sig, f, OK.1 , d) : sig in sigma];
        fInv := inverseAutomorphism(f);
        
        // Frobenius und Inverse auf OL fortsetzen
        f := continue_frobenius(f, OL, Prime(Zp));
        fInv := inverseAutomorphism(f, fInv);
        
        // Frobenius und Inverses auf OL2 fortsetzen,
        // dabei jeweils die Bilder f^(-i)(b) merken, wobei i=0..d-1 und b
        // Erzeuger der Erweiterung OL2/OL
        //vprint cohomTerm, 3: "Berechne Frobenius und sein Inverses";
        B := Reverse(generators_pad(OL2,OL));
        fB := [**];
        for b in B do
            f := continue_frobenius(f, Parent(b), Prime(Zp));
            fInv := inverseAutomorphism(f, fInv);
            seq := [b];
            for i in [1..d-1] do
                Append(~seq, fInv(seq[#seq]));
            end for;
            Append(~fB, seq);
        end for;
        
        //vprint cohomTerm, 3: "Erzeuge Fortsetzungen mit gewuenschtem Verhalten";
        // setze nun sigma forst
        // Falls sigma=f^j auf OK, j=0..d-1, und [OL2:OL]=m
        // Dann gilt fuer die Fortsetzungen sigmaHut|_Knr = F^(k*d+j)
        // Wir wollen die Fortsetzungen sigmaHut, so dass
        // F^(d-j)*sigmaHut|_Knr = id|_Knr
        // Wir waehlen also die Bilder zu Finv^(d-j)
        sigmaHut := [];
        for i in [1..#sigma] do
            sig := sigma[i];
            idx := frobIndex[i];
            // Liste f^(-idx)(b) fuer die Erzeuger von OL2/OL
            B := [* seq[((d-idx) mod d) + 1]   : seq in fB *];
            // sigma schrittweise fortsetzen
            for r in B do
                sig := map<Parent(r) -> Parent(r) | x:-> continued_automorphism_image(x, sig, r) >;
            end for;
            Append(~sigmaHut, sig);
        end for;
        
        // so konstruiert, dass das Inverse auf Knr einfach zu berechnen ist
        sigmaKnrExponentInv := frobIndex;
    end if;
    
    return sigmaHut, frobIndex, sigmaKnrExponentInv;
end intrinsic;

intrinsic inverseAutomorphism(aut::Map) -> Map
{ Compute the inverse of the given automorphism.
  If inv is given, the inverse is computed as a continuation of inv if possible. }
    local K;
    
    K := pAdicRing(Domain(aut));
    return inverseAutomorphism(aut, hom<K-> K| x:-> x>);
end intrinsic;

intrinsic inverseAutomorphism(aut::Map, inv::Map) -> Map
{ Compute the inverse of the given automorphism.
  If inv is given, the inverse is computed as a continuation of inv if possible. }
    local L, K, B, invAut, b, cont, i, c;
    
    L := Domain(aut);
    K := Domain(inv);
    
    invAut := inv;
    B := Reverse(generators_pad(L,K));
    for b in B do
        if (aut(b) eq b) then
            invAut := map<Parent(b) -> Parent(b) | x :-> continued_automorphism_image(x, invAut, b) >;
        else
            f := MinimalPolynomial(b);
            rts := {r[1] : r in Roots( Polynomial([ invAut(x) :  x in Coefficients(f)]), L) };
//"AliChanged";    
            //assert #rts eq Degree(f);
            cont := [ map< Parent(b) -> L | 
                x :-> &+([ invAut(elem[i])*r^(i-1) : i in [1..#elem]] 
                where elem is elem_to_seq(x,BaseRing(Parent(b)))) >  
                : r in rts 
            ];
            _, i := Maximum([Valuation(aut(c(b))-b) :  c in cont]);
            invAut := cont[i];
        end if;
    end for;
    
    return invAut;
end intrinsic;


intrinsic find_power(g, f, x, max)->.{}
    local p, m, fX, gX, X;
    
    if Type(x) eq SeqEnum then
        X := x;
    else
        X := [x];
    end if;
    
    
    m := (max eq 0 select Infinity() else max );
    p := 0;
    fX := X;
    gX := [g(x) : x in X];
    
    while p le m do
        if &and([ is_zero(gX[i]-fX[i]) : i in [1..#fX]]) then
            return p;
        end if;
        fX := [f(x) : x in fX];
        p := p+1;
    end while;
    
    error "Maximum power exceeded!!!";
end intrinsic;


/*
 * Given x=\sum a_i p^i, compute \sum m(a_i) b^i.
 */
intrinsic continued_automorphism_image(x, m, b)->.{}
    if Type(x) in {RngPadElt, FldPadElt} then
        if x eq ChangePrecision(Zero(Parent(x)),Precision(x)) then
            return ChangePrecision(Zero(Parent(x)),Precision(x));
        end if;
    end if;
    y := ElementToSequence(x);
    return &+([ m(y[i])*b^(i-1)  :  i in [1..#y] ]);
end intrinsic;


intrinsic elem_to_seq(x, K)->.{}
    if RingOfIntegers(Parent(x)) cmpeq RingOfIntegers(K) then
        return [x];
    end if;
    y := ElementToSequence(x);
    while not RingOfIntegers(Parent(y[1])) eq RingOfIntegers(K) do
        y := &cat([Coefficients(y[j]) : j in [1..#y]]);
    end while;
    return y;
end intrinsic;

intrinsic tupelProduct(t1, t2)->.{}
    return < t1[i]*t2[i]  : i in [1..#t1]>;
end intrinsic;

intrinsic tupelQuotient(t1, t2)->.{}
    return < t1[i]/t2[i]  : i in [1..#t1]>;
end intrinsic;

intrinsic tupelInverse(t)->.{}
    return < t[i]^(-1) : i in [1..#t] >;
end intrinsic;

intrinsic precompute_map(m)->.{}
    local domSeq, x, img;
    
    domSeq := [x : x in Domain(m)];
    img := [m(x) : x in domSeq];
    
    return map< Domain(m) -> Codomain(m) | x :->  img[Index(domSeq,x)] >;
end intrinsic;

/*
intrinsic maximal_unram_subext_simple(OL)->.{}
    local Zp, OK;
    
    if Type(OL) eq RngPad then
        Zp := pAdicRing(OL);
    else
        Zp := pAdicField(OL);
    end if;
    if isUnramified(OL,Zp) then
        OK := OL;
    elif isTotallyRamified(OL,Zp) then
        OK := Zp;
    else
        OK := OL;
        while AbsoluteDegree(OK) gt AbsoluteInertiaDegree(OK) do
            OK := BaseRing(OK);
        end while;
        assert AbsoluteDegree(OK) eq AbsoluteInertiaDegree(OK);
    end if;
    // habe nun Erweiterungen OL/OK/Zp 
    // mit OL/OK voll verzweigt und OK/Zp unverzweigt
    assert isTotallyRamified(OL, OK);
    assert isUnramified(OK, Zp);
    // OK soll "einfache" Erweiterung sein
    assert Degree(OK) eq Degree(OK, Zp);
    
    return OK;
end intrinsic;
*/
intrinsic isUnramified(L::FldPad, K::FldPad) -> BoolElt
{ Test whether L/K is unramified. }
    return InertiaDegree(L,K) eq Degree(L,K);
end intrinsic;
intrinsic isUnramified(L::RngPad, K::RngPad) -> BoolElt
{ Test whether L/K is unramified. }
    return InertiaDegree(L,K) eq Degree(L,K);
end intrinsic;

intrinsic isTotallyRamified(L::FldPad, K::FldPad) -> BoolElt
{ Test whether L/K is totally ramified. }
    return RamificationIndex(L,K) eq Degree(L,K);
end intrinsic;
intrinsic isTotallyRamified(L::RngPad, K::RngPad) -> BoolElt
{ Test whether L/K is totally ramified. }
    return RamificationIndex(L,K) eq Degree(L,K);
end intrinsic;

/* generates_maximal_unram(OL2::RngPad, OL::RngPad, OK::RngPad) -> BoolElt
 * Let OK be the maximal unramified extension of Zp in OL and OL2 an
 * extension of OL.
 * Test whether the maximal unramified extension OK2 in OL2/Zp can
 * be deduced from the extension OL2/OL, i.e. whether the coefficients of
 * the defining polynomial of OL2/OL are in OK.
 */
intrinsic generates_maximal_unram(OL2, OL, OK)->.{}
    local fields, E, d;
    
    if Degree(OL2, OL) eq 1 then
        return true;
    end if;
    
    //assert BaseRing(OL) eq OK;
    
    // erzeuge Koerperturm
    fields := Reverse(field_tower(OL2, OL));
    
    // OL selbst wollen wir nicht
    Remove(~fields,1);
    
    d := Degree(OL, OK);
    for E in fields do
        f := DefiningPolynomial(E);
        coeff := Coefficients(f);
        // Koeffizienten in Basisdarstellung
        coeff := [elem_to_seq(c, OK) : c in coeff];
        // jetzt darf jeweils nur jeder d-te Eintrag ungleich Null sein
        // entferne die Eintraege i*(d-1)+1
        for x in coeff do
            c := x;
            for i in [0..#c/d-1] do
                //print "entferne",  i*(d-1)+1;
                Remove(~c, i*(d-1)+1);
            end for;
            if SequenceToSet(c) ne {ChangePrecision(Zero(OK), Minimum([Precision(y) : y in c]))} then
                return false;
            end if;
        end for;
    end for;
    
    return true;
end intrinsic;

/* 
 * Compute the list of extensions from K to L, i.e. 
 * if L=E_1/E_2/.../E_n=K, return [*E_1,...,E_n*].
 */
intrinsic field_tower(L, K)->.{}
    lst := [* *];
    E := L;
    while not E eq K do
        Append(~lst, E);
        if E eq BaseRing(E) then
            error "K is not a subfield of L";
        end if;
        E := BaseRing(E);
    end while;
    Append(~lst, K);
    return lst;
end intrinsic;

intrinsic CLocalFundamentalClassSerre(L::FldPad, K::FldPad, steps::RngIntElt : psi:= 0, NoNorm := false) -> Map
{Compute the cocyle representing the local fundamental class of L/K
up to the given precision. Optionally, one can pass the map psi:G->Aut(L/K)
computed by AutomorphismGroup(L,K). If the parameter NoNorm=true,
the algorithm without norm equations is used.}
    return ali_cocycle_lfc_G(RingOfIntegers(L),RingOfIntegers(K),psi,steps,NoNorm);
end intrinsic;

intrinsic CLocalFundamentalClassSerre(L::RngPad, K::RngPad, steps::RngIntElt : psi:= 0, NoNorm := false) -> Map
{Compute the cocyle representing the local fundamental class of L/K
up to the given precision. Optionally, one can pass the map psi:G->Aut(L/K)
computed by AutomorphismGroup(L,K). If the parameter NoNorm=true,
the algorithm without norm equations is used.}
    return ali_cocycle_lfc_G(L,K,psi,steps,NoNorm);
end intrinsic;

intrinsic ClCocycleLFC(L::FldPad, K::FldPad, steps::RngIntElt : psi:= 0, NoNorm := false) -> Map
{Compute the cocyle representing the local fundamental class of L/K
up to the given precision. Optionally, one can pass the map psi:G->Aut(L/K)
computed by AutomorphismGroup(L,K). If the parameter NoNorm=true,
the algorithm without norm equations is used.}
    return ali_cocycle_lfc_G(RingOfIntegers(L),RingOfIntegers(K),psi,steps,NoNorm);
end intrinsic;

intrinsic ClCocycleLFC(L::RngPad, K::RngPad, steps::RngIntElt : psi:= 0, NoNorm := false) -> Map
{Compute the cocyle representing the local fundamental class of L/K
up to the given precision. Optionally, one can pass the map psi:G->Aut(L/K)
computed by AutomorphismGroup(L,K). If the parameter NoNorm=true,
the algorithm without norm equations is used.}
    return ali_cocycle_lfc_G(L,K,psi,steps,NoNorm);
end intrinsic;



/////////////////////CLocalBrauerGroup//////
 

 locBrGrp := recformat<
    L    : FldNum,        // the corresponding global field
    P    : RngOrdIdl,     // with prime ideal
    p    : RngIntElt,
    //m    : RngIntElt      // precision needed

    M    : GrpAb,        // the module
    actM : Map,          // group action on M
    qM   : Map,          // projection map L^\times ->> M
    theta: RngOrdElt,

    C    : ModCoho,      // local cohomology group structure
                         // computed for M and the subgroup H
    f1   : Map,          // corresponding map of M into the local structure
    lfc  : ModTupRngElt  // local fundamental class
>;

declare attributes FldNum: localBrauerGroups;


intrinsic AliLocalBrauerGroup(L::FldNum, p::RngIntElt : autMap := 0, lfc := false) -> Rec
{ Compute the local cohomology group at p: H^2(G_P,L_P^\times).
Returns a record containing:
  the cohomology structure C,
  the module M=L_P^\times/X s.t. H^2(G_P,L_P^\times)=H^2(G_P,M),
  the group action on M,
  a map qM:L->M,
  a map f1 from M onto the internal module,
  and the representant of the local fundamental class (if required).
}
    P := Factorization(p*RingOfIntegers(L))[1,1];
    return AliLocalBrauerGroup(L,P : autMap := autMap, lfc := lfc);
end intrinsic;



intrinsic AliLocalBrauerGroup(L::FldNum, P::RngOrdIdl : autMap := 0, lfc := false) -> Rec
{ Compute the local cohomology group at p: H^2(G_P,L_P^\times).
Returns a record containing:
  the cohomology structure C,
  the module M=L_P^\times/X s.t. H^2(G_P,L_P^\times)=H^2(G_P,M),
  the group action on M,
  a map qM:L->M,
  a map f1 from M onto the internal module,
  and the representant of the local fundamental class (if required).
}

    require IsAbsoluteField(L) and IsPrime(P) :
            "Absolute field L and prime ideal P required.";
    require IsNormal(L):
            "The field must be normal.";

    if not assigned L`localBrauerGroups then
        L`localBrauerGroups := [**];
    end if;
    p := Generator(P meet Integers());
    i := 0;
    for locBr in L`localBrauerGroups do
        i := i+1;
        if locBr`p cmpeq p then
            if not lfc or assigned locBr`lfc then
                return locBr;
            else
                // local fundamental class has to be computed
                // remove old record
                Remove(~L`localBrauerGroups, i);
            end if;
        end if;
    end for;

    if autMap cmpeq 0 then
        _,_,autMap := AutomorphismGroup(L);

  autMap := map< Domain(autMap) -> Codomain(autMap) | x :-> autMap(x^(-1)) >;
    end if;
    locBr := ali_local_brauer_group(L,P,autMap, lfc);
Append(~L`localBrauerGroups, locBr);
    return locBr;
end intrinsic;




intrinsic ali_local_brauer_group(L, P, psi, computeLFC)->.{}
    local G, pi, theta, m, LP, iota, H, psiL, lfc,
          V, M, proj, mm, HH, g, v, Y, mmY, y, X, qX, mmX, x,
          C, f1, f2, lfc2, H2;

    // compute prime ideal
    OL := RingOfIntegers(L);
    p := Generator(P meet Integers());
    vprint Brauer, 1: "Computing cohomology at", p;
    IndentPush();

    // Global Galois group
    G := Domain(psi);

    // Compute lattice
    pi := UniformizingElement(P);
    theta, m := lattice(P, pi, psi);
    vprint Brauer, 1: "lattice precision:" , m;
    vprint Brauer, 2: "theta: ", theta;

    // Localization with enough precision
    vprint Brauer, 1: "compute completion, prec=", 2*m+2;
    vtime Brauer, 1: LP, iota, psiL := completion_with_precision(L, P, psi, 2*m+2);
    // and Galois group
    //H, psiL := localized_automorphism_group(psi, P, iota, Automorphisms(LP, pAdicField(LP)));
    H := Domain(psiL);

    // compute V=LP^\times/X
   LIST:=[*LP,iota,psiL*];
   vprint Brauer, 1: "compute module";
    //X, mmX, qX := compute_LPmul_modX1(L, P, psi, LIST, theta, m);
    X, mmX, qX := compute_LPmul_modX(L, P, pi, psi, iota, LP, psiL, theta, m);
    // make this a right action
    mmX := map< H -> Aut(X) | g :-> mmX(g^(-1)) >;
    // map of L onto X
    qM := iota*qX;
    vprint Brauer, 2: "Dimension will be", #Invariants(X);

    // compute cohomology
    vprint Brauer, 1: "compute cohomology group";
    C := CohomologyModule(H, X, mmX);
    H2 := CohomologyGroup(C,2);
    // compute mappings to and from internal module V
    // f1: X -> V
    f1 := map< X -> RSpace(Integers(),Dimension(C)) |
        x:-> Eltseq(x),
        y:-> X!Eltseq(y)
    >;
// compute local fundamental class
    // with values in L^\times
    if computeLFC then
        vprint Brauer, 1: "compute cocycle, prec =",m;
       // lfc := cocycle_lfc_G(LP, pAdicField(LP), psiL, m, false);
     lfc := ali_cocycle_lfc_G(LP, pAdicField(LP), psiL, m, false);
        // read lfc in cohomology group
        vprint Brauer, 2: "identify cocycle";
        lfc2 := func< x | lfc(x[2]^(-1),x[1]^(-1)) @ qX @ f1 >;
        g := IdentifyTwoCocycle(C, lfc2);

        IndentPop();
        return rec< locBrGrp |
            L := L, P := P, p := p,
            M := X, actM := mmX, qM := qM,
            C:=C,f1:=f1, lfc := g,
            theta:=theta
        >;
    else
        IndentPop();
        return rec< locBrGrp |
            L := L, P := P, p := p,
            M := X, actM := mmX, qM := qM,
            C:=C,f1:=f1,
            theta:=theta
        >;
    end if;
end intrinsic;




////////////////////////////Added necessary files////////////////////////////////


intrinsic NormalBasisElement(OL :: RngOrd, h :: Map : rand := false) -> RngOrdElt
{ Given an order OL and a automorphism map h:G->Aut(OL),
  this function returns a normal basis element for OL. }
    local G, b, D, found;

    G := Domain(h);

    if not rand then
        for b in Basis(OL) do
            D := Matrix([ ElementToSequence( h(g)(b) ) : g in G ]);
            if Determinant(D) ne 0 then
                return OL!b;
            end if;
        end for;
    end if;

    found := false;
    while  not found  do
        b := OL ! [ Random(3) : i in [1..#G] ];
	D := Matrix([ ElementToSequence( h(g)(b) ) : g in G ]);
        if Determinant(D) ne 0 then
	   return OL!b;
	end if;
    end while;

    if not found then
        error "ERROR: No normal basis element found!!!";
    end if;

end intrinsic;


intrinsic maximal_unram_subext_simple(OL)->.{}
    local Zp, OK;

    if Type(OL) eq RngPad then
        Zp := pAdicRing(OL);
    else
        Zp := pAdicField(OL);
    end if;
    if isUnramified(OL,Zp) then
        OK := OL;
    elif isTotallyRamified(OL,Zp) then
        OK := Zp;
    else
        OK := OL;
        while AbsoluteDegree(OK) gt AbsoluteInertiaDegree(OK) do
            OK := BaseRing(OK);
        end while;
        assert AbsoluteDegree(OK) eq AbsoluteInertiaDegree(OK);
    end if;
    // habe nun Erweiterungen OL/OK/Zp 
    // mit OL/OK voll verzweigt und OK/Zp unverzweigt
    assert isTotallyRamified(OL, OK);
    assert isUnramified(OK, Zp);
    // OK soll "einfache" Erweiterung sein
    assert Degree(OK) eq Degree(OK, Zp);

    return OK;
end intrinsic;


intrinsic lattice(P, pi, psi)->.{}
    local OL;
    OL := Order(P);
    if Degree(OL) eq AbsoluteDegree(OL) then
        return lattice_absolute(P,pi,psi);
    else
        return lattice_relative(P,pi,psi);
    end if;
end intrinsic;


intrinsic lattice_check(P, psi)->.{}
    local OL;
    OL := Order(P);
    pi := UniformizingElement(P);
    if Degree(OL) eq AbsoluteDegree(OL) then
        return lattice_absolute(P,pi,psi);
    else
        return lattice_relative(P,pi,psi);
    end if;
end intrinsic;

intrinsic lattice_absolute(P, pi, psiL)->.{}
    local OL, p, theta, v, v1, erz, x, M, ZpGtheta, k, m, M1, M1I, j, M2, T;

    OL := Order(P);
    p := Generator(P meet Integers());
    G := Domain(psiL);
    rand := false;
//Ali changed here!!!

    repeat
        theta, v := lattice_generator_theta(psiL, P, pi, rand);
        rand := true;
        // erzeuger des Gitters global
        erz := [OL!(psiL)(g)(theta) : g in G];
        M := VerticalJoin( [Vector(ElementToSequence(x)) : x in erz ]);
        ZpGtheta := Lattice(M);
    until Rank(ZpGtheta) eq Degree(OL);

    // finde m mit P^m in ZpGtheta
    // einfacher Ansatz:
    k := Index(StandardLattice(Degree(OL)), ZpGtheta);
    m := Valuation(k*OL, P);
    //m := Valuation(k, setting`p)*RamificationIndex(setting`P);
    // "p=2 sometime its lengthy";
    if m eq 0 then
        return theta, v+1;

  end if;

    // kleinstes m
    // schreibe Basis von ZpGtheta in Matrix
    M1 := Matrix(Rationals(), [ElementToSequence(x) : x in erz]);
    M1I := M1^(-1);
    for j in [v+1..m] do
        m := j;
        // schreibe Basis von P^m in Matrix
        M2 := Matrix(Rationals(), [ElementToSequence(x) : x in Basis(P^j)]);
        // Basiswechsel durch T: M1*T=M2
        T := M1I*M2;
        // Elemente in T sollen nach Lokalisierung bei p ganz sein
        if not IsDivisibleBy(Denominator(T),p) then
            break;
        end if;
    end for;

    return theta, m;
end intrinsic;


intrinsic lattice_relative(P, pi, psiL)->.{}
    local OL,OK,p,G,theta,b,erz,g,M,x,y,ZpGtheta,k,m,M1,M1I,M2,j,T;

    OL := Order(P);
    OK := BaseRing(OL);
    assert(BaseRing(OK) cmpeq Integers());
    p := Generator(P meet Integers());
    G := Domain(psiL);
    rand := false;

    repeat
        theta, v := lattice_generator_theta(psiL, P, pi, rand);
        //print theta;
        rand := true;
        // erzeuger des Gitters global
        erz := [OL!(psiL)(g)(theta) : g in G];
        // multiplizieren mit Basis von OK
        erz := [x*(OL!y) : x in erz, y in Basis(OK)];
        M := VerticalJoin( [Vector(&cat([ Eltseq(y) : y in Eltseq(x)])) : x in erz ]);
        ZpGtheta := Lattice(M);
    until Rank(ZpGtheta) eq AbsoluteDegree(OL) and &and([x in Integers() : x in Eltseq(M) ]);

    // finde m mit P^m in ZpGtheta
    // einfacher Ansatz:
    //print M;
    k := Index(StandardLattice(AbsoluteDegree(OL)), ZpGtheta);
    m := Valuation(k*OL, P);

    if m eq 0 then
        return theta, v+1;
    end if;

    // kleinstes m
    // schreibe Basis von ZpGtheta in Matrix
    M1 := Matrix(Rationals(), M);
    M1I := M1^(-1);
    for j in [v+1..m] do
        m := j;
        // schreibe Basis von P^m in Matrix
        M2 := Matrix(Rationals(),  &cat([
            [ &cat([Eltseq(y) : y in Eltseq(b*x)]) : x in Basis(OK)]
            : b in Basis(P^j)
        ]));
        // Basiswechsel durch T: M1*T=M2
        T := M1I*M2;
        // Elemente in T sollen nach Lokalisierung bei p ganz sein
        if not IsDivisibleBy(Denominator(T),p) then
            break;
        end if;
    end for;

    return theta, m;
end intrinsic;


intrinsic lattice_generator_theta(psi, P, pi, rand)->.{}
    local OL, p, theta, v, v1;

    OL := Order(P);
    p := Generators(P meet Integers())[1];
   //checking here but doesn't work
  theta := NormalBasisElement(OL, psi : rand := rand);
   // theta := NormalBasisElement_check(OL, psi : rand := rand);
    v := Valuation(theta, P);
    v1 := 1+Floor(AbsoluteRamificationIndex(P)/(p-1));
    v := Maximum(0,v1-v);
    theta := OL!(theta*(pi)^v);

    return theta, v;
end intrinsic;

function completion_with_precision(L, P, psi, precision)
    local prec, min, err, compatible;

    if Generator(P meet Integers()) eq 2 then
        prec := Maximum(precision,50);
    else
        prec := Maximum(precision,50);
        //prec := Maximum(precision,30);
    end if;

    repeat
        err := false;
        //Ali changed
      compatible := false;
        try
            //print "compute completion", prec;
            LP, iota := Completion(L, P : Precision:=prec);
           // ChangePrecision(~LP,prec);
           // iota:=map<L->LP| x:-> LP!iota(x)>;
            autLP := Automorphisms(LP, pAdicField(LP));
            _, psiLP := localized_automorphism_group(psi, P, iota, autLP);
            /*H := [g : g in Domain(psi) | &and([  psi(g)(x) in P   : x in Generators(P)]) ];
            H := sub< Domain(psi) | H>;
            HH := [h : h in H];
            maps := [map< LP -> LP | x:-> x @@ iota @ psi(h) @ iota> :  h in HH];
            psiLP := map< H -> maps | h :-> maps[Index(HH,h)] >;
            */
            //print "test compatibility";
//"Ali changed to check for higher degree";

    min := Minimum(test_G_compatible_ali(iota, psi, psiLP, true, 0)
                join test_G_compatible_ali((iota)^(-1), psiLP, psi, false, P));
 /* min := Minimum(test_G_compatible(iota, psi, psiLP, true, 0)
                join test_G_compatible((iota)^(-1), psiLP, psi, false, P));
*/
              compatible := (min ge precision);

        catch e
            //print e`Object;
            err := true;
        end try;

        prec := 2*prec;
        if err then
            continue;
        end if;
    until (not err) and compatible;
//until (not err);
    return LP, iota, psiLP;
end function;



intrinsic completion_with_precision1(L, P, psi, precision)->.
{}
return completion_with_precision1(L, P, psi, precision);
end intrinsic;





intrinsic completion_with_precision(L, P, psi, precision)->.
{}
return completion_with_precision(L, P, psi, precision);
end intrinsic;



intrinsic localized_automorphism_group(m, P, iota, AutLoc)->.{}
    local G, GP, f, L, OL, Rts, RtsLok, i,j,prec,index, z, y, S;

    G := Domain(m);
    // Untergruppe von G
    //H := DecompositionGroup(P);
    GP := [g : g in G | &and([  m(g)(x) in P   : x in Generators(P)]) ];
    GP := sub< G | GP>;

    // Wenn G und H gleich sind, kann es sein, dass Magma die Gruppen
    // unterschiedlich aufzaehlt. D.h.
    // G eq H liefert true und
    // [g : g in G] eq [g : g in H] liefert false

    // dieses Verhalten ist nicht ganz nachvollziehbar und wird
    // hiermit umgangen
    if G eq GP then
        GP := G;
    end if;

    L := Domain(m(GP.1));
    OL := Domain(AutLoc[1]);
    Rts := Roots(DefiningPolynomial(L), L);
    //RtsLok := Roots(ChangePrecision(Polynomial(OL, ElementToSequence(f)),OL`DefaultPrecision));
    RtsLok := Roots(DefiningPolynomial(L), OL);

    assert #Rts eq #RtsLok;

    // Zuordnung:    globale Nst x <-> lokale Nst y
    z := [];
    for i in [1..#Rts] do
        S := [ Valuation(RtsLok[j,1] - OL!iota(Rts[i,1])) : j in [1..#RtsLok] ];
        prec, index := Maximum(S);
        z := z cat [index];
    end for;
    //print z;

    // Zuordnung:    g in AutLoc <-> index von g(RtsLok[1]) in RtsLok
    y := [];
    for i in [1..#AutLoc] do
        S := [ Valuation(AutLoc[i](RtsLok[1,1]) - RtsLok[j,1] ) : j in [1..#RtsLok] ];
        //print S;
        prec, index := Maximum(S);
        y := y cat [index];
    end for;
    //print y;

    // Zuordnung:    Index globale Nst x  <->  Index von g in AutLoc, so dass g(RtsLok[1])=y
    z := [ Index(y, z[i]) : i in [1..#z] ];

    return GP, map< GP -> AutLoc | x :-> local_map(x, m, AutLoc, Rts, z) >;
end intrinsic;













intrinsic test_G_compatible_ali(phi, actDom, actCodom, modP, prime)->.{}
    local D, B, gens, actD, actB, seq, U;

    if Type(prime) eq RngOrdIdl then
        modP := true;
    end if;

    D := Domain(phi);
    B := Codomain(phi);

    if Type(Domain(actDom)) ne SetCart then
        G := Domain(actDom);
        actD := map< car<G, D> -> D | x :-> actDom(x[1])(x[2]) >;
    else
        G := Component(Domain(actDom),1);
        actD := actDom;
    end if;

    if Type(Domain(actCodom)) ne SetCart then
        H := Domain(actCodom);
        //assert G eq Domain(actCodom);
        actB := map< car<H, B> -> B | x :-> actCodom(x[1])(x[2]) >;
    else
        H := Component(Domain(actCodom),1);
        //assert G eq Component(Domain(actCodom),1);
        actB := actCodom;
    end if;

    if G eq H then
        // groups equal
        U := G;
    else
        // take the smaller group
        if #H lt #G then
            U := H;
        else
            U := G;
        end if;
        // and make sure the elements can be read in the other group
       // assert &and([x in G and x in H :x in U]);
    end if;


    if Type(D) in {RngOrd, FldNum, ModTupRng, ModTupFld} then
        gens := Basis(D);
    elif Type(D) in {FldPad, RngPad} then
       gens := basis_generators(D);
      // gens := AbsoluteBasis(D);
    else
        print "not yet implemented: Generators/Basis for ", Type(D);
        try
            gens := Basis(D);
        catch e
            gens := Generators(D);
        end try;
    end if;
   //Ali changed for faster 
    if modP then
    seq:=[];
    for x in gens do
       for sig in Generators(U) do
            Append(~seq,phi(actD(sig, x)) - actB(sig, phi(x)));
       end for;
    end for;
   //  seq := [ phi(actD(sig, x)) - actB(sig, phi(x)) : x in gens, sig in U];
        if Type(B) in {FldPad,RngPad} then
            return {Valuation(x) : x in seq};
        elif Type(B) in {FldNum,RngOrd} then
            if Type(prime) ne RngOrdIdl then
                error "Prime Ideal for Valuation needed!";
            end if;
            return {Valuation(x, prime) : x in seq};
        else
            error "not yet implemented: Valuation";
        end if;
    else
        //seq := [ [* sig, x, B!(phi(actD(sig, x)) - actB(sig, phi(x))) *] : x in gens, sig in G ];
        //print seq;
        return &and({ phi(actD(sig, x)) eq actB(sig, phi(x)) : x in gens, sig in U});
    end if;
end intrinsic;


intrinsic test_G_compatible(phi, actDom, actCodom, modP, prime)->.{}
    local D, B, gens, actD, actB, seq, U;

    if Type(prime) eq RngOrdIdl then
        modP := true;
    end if;

    D := Domain(phi);
    B := Codomain(phi);

    if Type(Domain(actDom)) ne SetCart then
        G := Domain(actDom);
        actD := map< car<G, D> -> D | x :-> actDom(x[1])(x[2]) >;
    else
        G := Component(Domain(actDom),1);
        actD := actDom;
    end if;

    if Type(Domain(actCodom)) ne SetCart then
        H := Domain(actCodom);
        //assert G eq Domain(actCodom);
        actB := map< car<H, B> -> B | x :-> actCodom(x[1])(x[2]) >;
    else
        H := Component(Domain(actCodom),1);
        //assert G eq Component(Domain(actCodom),1);
        actB := actCodom;
    end if;

    if G eq H then
        // groups equal
        U := G;
    else
        // take the smaller group
        if #H lt #G then
            U := H;
        else
            U := G;
        end if;
        // and make sure the elements can be read in the other group
       // assert &and([x in G and x in H :x in U]);
    end if;
 if Type(D) in {RngOrd, FldNum, ModTupRng, ModTupFld} then
        gens := Basis(D);
    elif Type(D) in {FldPad, RngPad} then
      // gens := basis_generators(D);
       gens := AbsoluteBasis(D);
    else
        print "not yet implemented: Generators/Basis for ", Type(D);
        try
            gens := Basis(D);
        catch e
            gens := Generators(D);
        end try;
    end if;
   //Ali changed for faster
    if modP then
    seq:=[];
    for x in gens do
       for sig in U do
            Append(~seq,phi(actD(sig, x)) - actB(sig, phi(x)));
       end for;
    end for;
   //  seq := [ phi(actD(sig, x)) - actB(sig, phi(x)) : x in gens, sig in U];
        if Type(B) in {FldPad,RngPad} then
            return {Valuation(x) : x in seq};
        elif Type(B) in {FldNum,RngOrd} then
            if Type(prime) ne RngOrdIdl then
                error "Prime Ideal for Valuation needed!";
            end if;
            return {Valuation(x, prime) : x in seq};
        else
            error "not yet implemented: Valuation";
        end if;
    else
        //seq := [ [* sig, x, B!(phi(actD(sig, x)) - actB(sig, phi(x))) *] : x in gens, sig in G ];
        //print seq;
        return &and({ phi(actD(sig, x)) eq actB(sig, phi(x)) : x in gens, sig in U});
    end if;
end intrinsic;





/*
intrinsic basis_generators(L::FldPad)-> SeqEnum
    {it returns only the uniforisers of L and its base fields}
  gen :=[];
   repeat
       Append(~gen,L.1);
       l:= BaseField(L);
       if Degree(l,PrimeField(l)) gt 1 then
	  Append(~gen,L!l.1);
       end if;	  
   until Degree(l, PrimeField(l)) eq 1;
return gen;
end intrinsic;*/

intrinsic basis_generators(L::FldPad)-> SeqEnum
    {it returns only the uniforisers of L and its base fields}
  gen :=[];
   repeat
       Append(~gen,L.1);
       L:= BaseField(L);
       /*if Degree(l,PrimeField(l)) gt 1 then
          Append(~gen,L!l.1);
       end if;  */
   until Degree(L, PrimeField(L)) eq 1;
return gen;
end intrinsic;

intrinsic basis_generators(L::RngPad)-> SeqEnum
   {it returns only the uniforisers of L and its base fields}
  gen :=[];
   repeat

       Append(~gen,L.1);
       L:= BaseRing(L);
   /*    if Degree(l,PrimeRing(l)) gt 1 then
          Append(~gen,L!l.1);
       end if;*/
   until Degree(L, PrimeRing(L)) eq 1;
return gen;
end intrinsic;


intrinsic local_map(g, m, HomL, Rts, z)->.{}
//localMap(g::., m::Map, HomL::SeqEnum, Rts::SeqEnum, z::SeqEnum) -> Map
    // Nehme die globale Nst x0, die auf die erste lokale Nst abbildet
    first := Index(z,1);
    x := m(g)(Rts[first,1]);
    // Finde Index der globalen Nst x, so dass g(x0)=x
    S := [ x- Rts[i,1] : i in [1..#Rts] ];
    j := Index(S, 0);
    // Der Index der lokalen Abb, die das gleiche tut, steht in z
    return HomL[z[j]];
end intrinsic;

