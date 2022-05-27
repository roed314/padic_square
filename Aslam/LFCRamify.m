
////////////////////////make this for tower of extensions////////////////////////////////////////////////////////////////////////////////////////////////////

declare verbose Frob_Eq, 2;


///////////////////Frobenius Equation uniding Hilbert 90 theorem////////////////////////////


intrinsic solve_frob_linear(L,A,n)->.{}
  phi :=FrobeniusAutomorphism(L);
  prec := n;
  U,mU := UnitGroup(L:Prec:= prec);
  gens := CUnitGroupGenerators(L);
  F := FreeAbelianGroup(#gens);
  h := hom<F->U| [(phi(x) div x)@@mU : x in gens]>;
   A_seq:=[ Eltseq(a@@mU@@h): a in A];
  B := [];
  for j in [1..#A] do
     Append(~B, &*[gens[i]^A_seq[j,i]: i in [1..#gens]] );
//  B := &*[gens[i]^A_seq[i]: i in [1..#gens]];
   end for;
return B;
end intrinsic;

 
 intrinsic Hilbert_Coeff(c,phi)->.{}
         L := Domain(phi);
         pi := UniformizingElement(L);
         basis := Basis(L) ;
         gamma := basis[2];
	 d := InertiaDegree(L);
         B := [];
	 for i in [0..d-1] do
             Append(~B,gamma@(phi^(i)));
         end for;
return B;
end intrinsic;



/*
 * Creates a list of generators of L/K.
 * If L=L_n/L_{n-1}/.../L_0=K, this returns the list [L_n.1,L_{n-1}.1,...,L_0.1].
 */
function generators_pad(L, K)
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
end function;
 
 /*
intrinsic frobenius_equation_Hilbert_ppower(set,phi,B)->.{time consuming bcz of operations of Frobenius automorphisms}

	//B := [];
	L := Domain(phi);
	p := Prime(L);
        pi := UniformizingElement(L);
        basis := Basis(L) ;
        gamma := B[1];
	 
	 d := InertiaDegree(L);
         quot := Floor(d div p);
         rem := d-p*quot;  
	c := set[1]^-1;        
 //        A := [gamma,c*gamma@phi];
          aa:=[c,c@phi];// c^{1+phi+..+ phi^(p-1)}
          for j in [3..p] do 
	      Append(~aa,aa[j-1]@phi);
	   end for;
 	  cp := &*( aa );
 	  A :=[gamma, c*(gamma@phi)];
	    
	  for j in [3..quot] do
	    Append(~A, c*(A[j-1]@phi))  ;
           end for;
          term := &+(A);
          B := [term];
          CP := [cp@phi^(i*p): i in [0..quot-1]];
 	for j in [2..quot] do
            Append(~B,&*(CP[1..j-1])*(B[j-1]@( phi^p )) );
          end for; 
 
            for j in [3..p] do
            Append(~A,c* A[j-1]@phi);
            end for;
	 sum := &+(A);
	for i in [2..quot] do 
            Append(~A, )

// for i in [0..d-1] do
	//     Append(~B,gamma@(phi^(i))); 
	 //end for; 
	seq := [];
	for i in [1..#set] do
	    c := set[i]^-1;
	    A := [1,c];
	    for i in [3..d] do
	    Append(~A,A[i-1]*c@(phi^(i-2))) ;
	    end for;
	    b := &+([A[i]*B[i]: i in [1..d]]);
	   assert Norm(b) ne 0 ;
	   Append(~seq, b);
	   end for;
return seq;
end intrinsic;
	   */ 
intrinsic frobenius_equation_Hilbert_old(set,phi)->.{time consuming bcz of operations of Frobenius automorphisms}
 L := Domain(phi);
        pi := UniformizingElement(L);
        basis := Basis(L) ;
        gamma := 1+basis[2];
        d := InertiaDegree(L);
    seq := [];
    for i in [1..#set] do
        c := set[i]^-1;
        A := [gamma,c*gamma@phi];
        for j in [3..d] do
            Append(~A,c* A[j-1]@phi);
        end for;
        c_frob := &+(A);
        if Norm(c_frob) eq 0 then
           c_frob := Hilbert_Ninety1(c,phi);
         end if;
        Append(~seq, c_frob);
    end for;
return seq;
end intrinsic;



intrinsic frobenius_equation_Hilbert(set,phi)->.{This solves the frobenius equations and is good enought because of iteration process}
 L := Domain(phi);
        pi := UniformizingElement(L);
        basis := Basis(L) ;
        gamma := 1+basis[2];
	d := InertiaDegree(L);
    seq := [];
    for i in [1..#set] do
	c := set[i]^-1;
	A := [gamma,c*gamma@phi];
         //"d is the order of phi which is Inertia degree pf OL and phi =phi'^{d'}";
	for j in [3..d] do
	    Append(~A,c* A[j-1]@phi);
	end for;
	c_frob := &+(A);
	if Norm(c_frob) eq 0 or Valuation(c_frob) ge 1 then 
	   c_frob := Hilbert_Ninety1(c,phi,basis);
	 end if;  
	Append(~seq, c_frob);
    end for;	
return seq;
end intrinsic;

intrinsic Hilbert_Ninety1(c,phi,basis)->.{Computes an element x satisfying x^(phi -1)=c using Hilbert 90 formula twice at least beta neq 0}


        L := Domain(phi);
        pi := UniformizingElement(L);
        if Valuation(phi(pi)/pi - c) ge Precision(c)  then
           return pi;
        end if;
        d := InertiaDegree(L);
     //	c := c^(-1); 
       // basis := Basis(L) ;
        gamma := 1+basis[#basis];
        A:= [gamma,c*gamma@phi];
	for i in [3..d] do
            Append(~A,c* A[i-1]@phi);
        end for;

        c_frob := &+(A);
	
      if Norm(c_frob) eq 0 or Valuation(c_frob) ge 1 then 
	  gamma := Random(Integers(L));
          A:= [gamma,c*gamma@phi];
           for i in [3..d] do
              Append(~A,c* A[i-1]@phi);
          end for;
          return &+(A);
      else return &+(A);
      end if;
end intrinsic;	
 





intrinsic Hilbert_Ninety_old(c,phi)->.{Computes an element x satisfying x^(phi -1)=c}
        
 
	L := Domain(phi);
	pi := UniformizingElement(L);
        if Valuation(phi(pi)/pi - c) ge Precision(c)  then 
	   return pi;
	end if;   
	d := InertiaDegree(L);
        basis := Basis(L) ;
        gamma := basis[2];
        B:= [gamma,c*gamma@phi];
        A := [1,c];
        for i in [3..d] do
            Append(~A, A[i-1]*c@(phi^(i-2)));
            Append(~B, A[i]*gamma@(phi^(i-1)));
        end for;
       if &+(B) eq 0 then
          gamma := Random(L);
          B:= [1,c*gamma];
          A := [1,c];
          for i in [3..d] do
            Append(~A, A[i-1]*c@(phi^(i-2)));
            Append(~B, A[i]*gamma@(phi^(i-1)));
        end for;
        assert &+(B) ne 0;
        return &+B;
 else   return &+B;
 end if;
end intrinsic;



intrinsic LFC_Using_Hilbert(L,K,precision)->.{canonical class using fast norm equation and solving the Frobenius equation usin Hilbert 90}

    _,psi,_ := AutomorphismGroup(L,K);
        psi := map<Domain(psi) -> Codomain(psi) | g :-> psi(g^(-1))>;
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
if steps gt L`DefaultPrecision then
        error "Precision of p-adic field L not high enough to compute the cocycle!";
    end if;
    if steps gt Minimum([ Precision(r[1]) : r in Roots(DefiningPolynomial(L),L)]) then
        error "Precision of automorphisms of L not high enough to compute the cocyle!";
    end if;

          L1 := L;
          K1 := K;
	  OK := RingOfIntegers(K);
          e := RamificationIndex(L1,K1);
          OL1 := RingOfIntegers(L1);
          OL2 := ext<OL1 |e>;
	  u := UniformizingElement(K1)/UniformizingElement(L1)^e;
          vtime  Frob_Eq, 1: gamma1 := MyNormEquation_prec(OL2,OL1!u);
	 // vtime   CocycleLFC, 1: gamma1 := ClNormEquation(OL2,OL1!u);
	  pi:= gamma1*UniformizingElement(L1);
	  pi := OL2!pi;
//	  phi := FrobeniusAutomorphism(OL2, OL);//"this is expensive in large degree so find to accelerate";
          phi :=FrobeniusAutomorphism(OL2);
	  d :=InertiaDegree(L,K);
         if AbsoluteDegree(K) eq 1 then
             frobAction, GAction, frobIndex := galois_act_L_tensor_Knr(OL, OL2, psi, phi);
else
           frobAction, GAction, frobIndex := galois_act_L_tensor_Knr_ramify(OL,OK, OL2, psi, phi);
    end if;
      vtime Frob_Eq, 2:  pi_sigma := [GAction(g,<pi : i in [1..d]>)[1] : g in G];
          pisigmapi := [ OL2!(pi_sigma[i]/pi) : i in [1..#pi_sigma]];
          vprintf CocycleLFC, 1: "Solve Frobenius equations... ";
         // B := Hilbert_Coeff(pisigmapi[2],phi);
          u_sigma := [UniformizingElement(OL2)]; 
 vtime Frob_Eq, 2:  u_sigma1 := frobenius_equation_Hilbert(pisigmapi[2..#pisigmapi],phi);
	 //u_sigma1 := frobenius_equation_Hilbert(pisigmapi[2..#pisigmapi],phi);
         for i in [1..#u_sigma1] do 
	     Append(~u_sigma, u_sigma1[i]);
	 end for;    
/*	 for i in [1..#pisigmapi] do 
	     Append(~u_sigma, Hilbert_Ninety(pisigmapi[i]^(-1),phi));
	 end for; 
 */
	 //vtime  CocycleLFC, 1: u_sigma, phi := FrobeniusEquation(pisigmapi, steps, OL);
        if GetVerbose("CocycleLFC") ge 2 then
           vprint CocycleLFC, 2: "Test FrobeniusEquation result";
           assert &and({Valuation(phi(u_sigma[i])/u_sigma[i] - pisigmapi[i]) ge steps : i in [1..#u_sigma]});
         end if;
    // Kozykel
//    d := InertiaDegree(OL,Zp);
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

//    c := precompute_map(c);
     G_cart := [<x,y>: x,y in Generators(G)];
      //assert Minimum([ Minimum([Valuation(y[1]-y[i]) : i in [1..#y]])
        //where y is c(x) :  x in Domain(c)]) ge (precision+1);
   assert Minimum([ Minimum([Valuation(y[1]-y[i]) : i in [1..#y]]) where y is c(x) :  x in G_cart]) ge (precision+1);
 /*    if Degree(Codomain(c)[1]) ge 2 then
        // erste Komponente in L modulo pi^(precision+1)
        assert Minimum([ Minimum([ Valuation(z) : z in ElementToSequence(y[1])[2..Degree(Parent(y[1]))]])
            where y is c(x)  :  x in Domain(c)]) ge (precision+1);
    end if;*/
if Degree(Codomain(c)[1]) ge 2 then
        // erste Komponente in L modulo pi^(precision+1)
        assert Minimum([ Minimum([ Valuation(z) : z in ElementToSequence(y[1])[2..Degree(Parent(y[1]))]])
            where y is c(x)  :  x in G_cart]) ge (precision+1);
    end if;
    gamma := map< Domain(c) -> FieldOfFractions(L) | x :->  ( elem_to_seq(c(x)[1], L)[1] )^(-1) >;

    return gamma, psi;
end intrinsic;



 intrinsic precompute_map_iterate(m)->.{}
    local domSeq, x, img;

    domSeq := [x : x in Domain(m)];
    img := [];
    for x in domSeq do
        Append(~img, x@m);
    end for;
   // img := [m(x) : x in domSeq];

    return map< Domain(m) -> Codomain(m) | x :->  img[Index(domSeq,x)] >;
end intrinsic;







intrinsic maximal_unram_subext_ramify(OL,base)->.
{}    local Zp, OK;

   /* if Type(OL) eq RngPad then
        Zp := pAdicRing(OL);
    else
        Zp := pAdicField(OL);
    end if;*/

    if isUnramified(OL,base) then
        OK := OL;
    elif isTotallyRamified(OL,base) then
        OK := base;
    else
        OK := OL;
        while Degree(OK,base) gt InertiaDegree(OK,base) do
            OK := BaseRing(OK);
        end while;
        assert Degree(OK,base) eq InertiaDegree(OK,base);
    end if;
    // habe nun Erweiterungen OL/OK/Zp
    // mit OL/OK voll verzweigt und OK/Zp unverzweigt
    assert isTotallyRamified(OL, OK);
    assert isUnramified(OK, base);
    // OK soll "einfache" Erweiterung sein
    if isTotallyRamified(OL,BaseRing(OL)) and isTotallyRamified(BaseRing(OL),base) then
        assert Degree(OK) eq Degree(base);
end if;
 //assert Degree(OK) eq Degree(OK, base);

    return OK;
end intrinsic;



intrinsic galois_act_L_tensor_Knr_ramify(OL,base, OL2, psi, phi)->.{alternative of galois_act_L_tensor_Knr}
        if isTotallyRamified(OL,BaseRing(OL)) and isTotallyRamified(BaseRing(OL),base) then
           OK := maximal_unram_subext_ramify(OL,base);
         else
            OK := maximal_unram_subext_simple(OL);
         end if;
         Zp := base;
     if not generates_maximal_unram(OL2, OL, OK) then
        error "The maximal unramified extension in OL2/Zp cannot be deduced from OL2/OL!";
    end if;

    // Initialization
    G := Domain(psi);
    GG := [g : g in G];
    d := InertiaDegree(OL,Zp);

    // Compute i such that sigma_OK = phi^i for all sigma in G
    // and extensions sigmaHut of sigma such that sigmaHut^(-1)=phi^(-i) on K2
    sigmaHut, frobIndex := continuations_with_unram_restriction_ramify(OL,base, G, psi, OL2);

    // Frobenius automorphism on \prod L2
    L2 := FieldOfFractions(OL2);
    prodL2 := CartesianProduct([L2 : y in [1..d] ]);
    frobeniusMap := map< prodL2 -> prodL2 | x :->  < i eq d select phi(x[1]) else x[i+1]  : i in [1..d]> >;

    // action of G on \prod L2
    Gaction := map< car<G, prodL2> -> prodL2 | x :->
        apply_map( (frobeniusMap^((d-frobIndex[Index(GG,x[1])]) mod d))( x[2] ), sigmaHut[Index(GG, x[1])]) >;

    return frobeniusMap, Gaction, frobIndex;
end intrinsic;
intrinsic LFC_ramify_tower(L, K, precision)->.{computes in many tower subextensions}
 local steps, G, OL, Zp, Qp, pi, pi_sigma, pisigmapi, g, i, u_sigma, phi,
          e, OL2ub, OL2, OL1, L1, K1, incl, GG, AutL, AutL1, sigma, psi1,
          OLmal, m, u, bool, gamma1;
//if psi cmpeq 0 then
        _,psi,_ := AutomorphismGroup(L,K);
        psi := map<Domain(psi) -> Codomain(psi) | g :-> psi(g^(-1))>;
  //  end if;
is_equal := func< x,y | ChangePrecision(x,m) eq ChangePrecision(y,m)
                        where m := Minimum(Precision(x), Precision(y)) >;

steps := precision+2;
    G := Domain(psi);
    OL := RingOfIntegers(L);
    OK := RingOfIntegers(K);
    Zp := pAdicRing(RingOfIntegers(OL));
    Qp := FieldOfFractions(Zp);
    d := InertiaDegree(OL,Zp);

if d eq Degree(L,K) then
   //gamma := CLocalFundamentalClassSerre_check(L,K,Precision);
   gamma := CLocalFundamentalClassSerre(L,K,Precision);
end if;
      d := InertiaDegree(OL,OK);
      L1 := L;
      K1 := K;
      e := RamificationIndex(L1,K1);
      OL1 := RingOfIntegers(L1);
      OL2 := ext<OL1 |e>;
      u := UniformizingElement(K1)/UniformizingElement(L1)^e;
      vprintf CocycleLFC, 1: "Solve Norm equation... ";
       // vtime   CocycleLFC, 1: bool, gamma1 := MyNormEquation(OL2,m,OL1!u);
      vtime   CocycleLFC, 1: gamma1 := ClNormEquation(OL2,OL1!u);
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
        frobAction, GAction, frobIndex := galois_act_L_tensor_Knr_ramify(OL,OK, OL2, psi, phi);
        pi_sigma := [GAction(g,<pi : i in [1..d]>)[1] : g in G];
        pisigmapi := [ OL2!(pi_sigma[i]/pi) : i in [1..#pi_sigma]];
        vprintf CocycleLFC, 1: "Solve Frobenius equations... ";
        vtime  CocycleLFC, 1: u_sigma, phi := FrobeniusEquation(pisigmapi, steps, OL);


 if GetVerbose("CocycleLFC") ge 2 then
        vprint CocycleLFC, 2: "Test FrobeniusEquation result";
        assert &and({Valuation(phi(u_sigma[i])/u_sigma[i] - pisigmapi[i]) ge steps : i in [1..#u_sigma]});
    end if;

    // Kozykel
   // d := InertiaDegree(OL,Zp);
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
   assert Minimum([ Minimum([Valuation(y[1]-y[i]) : i in [1..#y]])
        where y is c(x) :  x in Domain(c)]) ge (precision+1);
    if Degree(Codomain(c)[1]) ge 2 then
        // erste Komponente in L modulo pi^(precision+1)
        assert Minimum([ Minimum([ Valuation(z) : z in ElementToSequence(y[1])[2..Degree(Parent(y[1]))]])
            where y is c(x)  :  x in Domain(c)]) ge (precision+1);
    end if;
 gamma := map< Domain(c) -> FieldOfFractions(L) | x :->  ( elem_to_seq(c(x)[1], L)[1] )^(-1) >;

    return gamma,psi;
end intrinsic;





intrinsic continuations_with_unram_restriction_ramify(OL,base, G, psi_OL_Zp, OL2)->.{this computes over unramified extension}
local Hom_OL_Zp, OK, Zp, sigma, sig, f,
          sigmaHut, frobIndex, sigmaKnrExponentInv;
Hom_OL_Zp := Codomain(psi_OL_Zp);
    OK := maximal_unram_subext_ramify(OL,base);
    Zp := PrimeRing(OK);
 if not generates_maximal_unram(OL2, OL, OK) then
        print "ACHTUNG: Fuer den angegebenen Koerperturm ist der schnelle Algorithmus nicht anwendbar!";
        print "Benutze den langsamen Algorithmus via primitiven Elementen.";

        error "Could not compute 'unramified behaviour' of extension";
    end if;
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

        f := FrobeniusAutomorphism(OK, Zp);
        frobIndex := [find_power(sig, f, OK.1 , InertiaDegree(OK, pAdicRing(OK))) : sig in sigma];
        fInv := inverseAutomorphism(f);

        sigmaHut := sigma;
        sigmaKnrExponentInv := [(d-i) mod d : i in frobIndex];

 else
        // Verhalten von sigma\in\Aut(OL,Zp) auf OK bestimmen
        f := FrobeniusAutomorphism(OK, Zp);
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

intrinsic maximal_unram_subext_simple_ramify(OL,Zp)->.{}
    //local Zp, OK;

   /* if Type(OL) eq RngPad then
        Zp := pAdicRing(OL);
    else
        Zp := pAdicField(OL);
    end if;*/
    if isUnramified(OL,Zp) then
        OK := OL;
    elif isTotallyRamified(OL,Zp) then
        OK := Zp;
    else
        OK := OL;
        while Degree(OK,Zp) gt InertiaDegree(OK,Zp) do
            OK := BaseRing(OK);
        end while;
     // assert AbsoluteDegree(OK) eq AbsoluteInertiaDegree(OK);
       assert Degree(OK,Zp) eq InertiaDegree(OK,Zp);
     end if;
    // habe nun Erweiterungen OL/OK/Zp 
    // mit OL/OK voll verzweigt und OK/Zp unverzweigt
    assert isTotallyRamified(OL, OK);
    assert isUnramified(OK, Zp);
    // OK soll "einfache" Erweiterung sein
   // assert Degree(OK) eq Degree(OK, Zp);

    return OK;
end intrinsic;


intrinsic galois_act_L_tensor_Knr_ramify(OL,base, OL2, psi, phi)->.{}
    local OK, Zp, G, GG, g, d, sigmaHut, frobIndex, L2, prodL2,
        frobeniusMap, Gaction, x;


    // Test whether algorithm is applicable
    OK := maximal_unram_subext_simple_ramify(OL,base);
  //  OK := BaseRing(OL);
    Zp := BaseRing(OK);

    if not generates_maximal_unram(OL2, OL, OK) then
        error "The maximal unramified extension in OL2/Zp cannot be deduced from OL2/OL!";
    end if;

    // Initialization
    G := Domain(psi);
    GG := [g : g in G];
   // d := InertiaDegree(OL,Zp);//working only over K;
     d := InertiaDegree(OL,OK);
    // Compute i such that sigma_OK = phi^i for all sigma in G
    // and extensions sigmaHut of sigma such that sigmaHut^(-1)=phi^(-i) on K2
    sigmaHut, frobIndex := continuations_with_unram_restriction_ramify(OL,base, G, psi, OL2);

    // Frobenius automorphism on \prod L2
    L2 := FieldOfFractions(OL2);
    prodL2 := CartesianProduct([L2 : y in [1..d] ]);
    frobeniusMap := map< prodL2 -> prodL2 | x :->  < i eq d select phi(x[1]) else x[i+1]  : i in [1..d]> >;

    // action of G on \prod L2
    Gaction := map< car<G, prodL2> -> prodL2 | x :->
        apply_map( (frobeniusMap^((d-frobIndex[Index(GG,x[1])]) mod d))( x[2] ), sigmaHut[Index(GG, x[1])]) >;

    return frobeniusMap, Gaction, frobIndex;
end intrinsic;

intrinsic continuations_with_unram_restriction_ramify(OL,base, G, psi_OL_Zp, OL2)->.{}
    local Hom_OL_Zp, OK, Zp, sigma, sig, f,
          sigmaHut, frobIndex, sigmaKnrExponentInv;

    //require Type(Domain(psi_OL_Zp(G.1))) in {RngPad, FldPad} :
    //        "Bad argument types\nAutomorphisms of p-adic rings/fields needed.";

    Hom_OL_Zp := Codomain(psi_OL_Zp);
    OK := maximal_unram_subext_simple_ramify(OL,base); //"already computed in earlier";
   // OK := base;
   // Zp := BaseRing(OK);//changed;
      Zp := PrimeRing(OK);
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






