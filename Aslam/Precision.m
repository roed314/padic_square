import "brauer.m"  : lattice_generator_theta, lattice, find_included_P_power_absolute;
     import "brauer.m"  :   compute_LPmul_modX;

          import "AliLFC.m"  : ali_cocycle_lfc_G;
  //  import "RelativeExtension.m"  :   CLocalFundamentalClassSerre_check;

GFCcomp := recformat<
    CohL : ModCoho,         // cohomology group
    f1CL : Map,             //   with map
    gfcId : ModTupRngElt,   // global fundamental class
    CL : GrpAb,             // module C_L
    psiCL : Map,            //   with G-action
    qCL : Map,              // projection J_L-->>C_L
    primes : SeqEnum,       // set of primes
    US : GrpAb,             // S-Units w.r.t primes
    mUS : Map,              //   with map to L
    kappaInf : SeqEnum,     // inclusions US --> ind^G US
    RepInf : SeqEnum,       //   with corresponding system of representatives
    inclJL : SeqEnum,       // inclusions J_{L_v} --> J_L
    inclUSJL : Map,
    lat : Any,              // lattice
    theta : SeqEnum        // lattice generators
>;



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



intrinsic frobenius_Hilbert_matrix(set,phi)->.{This solves the frobenius equations and is good enought because of iteration process}
 L := Domain(phi);
        pi := UniformizingElement(L);
        basis := Basis(L) ;
        gamma := 1+basis[2];
        //gamma :=1+ UniformizingElement(L);
      // gamma := ChangePrecision(gamma,60); 
       d := InertiaDegree(L);
   /* seq := [];
    for i in [1..#set] do
        c := set[i]^-1;
        A := [gamma,c*gamma@phi];*/
        D := [gamma@phi];
    for j in [2..d] do
        Append(~D,D[j-1]@phi);
    end for;
     D:=[L!1: i in [1..d]];
     D:= Transpose(Matrix([D]));
     S:=[[1,set[i]^-1]: i in [1..#set]];
   //  S:=[[1,ChangePrecision(set[i]^-1,60)]: i in [1..#set]];
     for j in [1..d-2] do
         for n in [1..#S] do
             Append(~S[n], S[n,2]*S[n,j+1]@phi);
          end for;
     end for;
   A :=[Matrix([a]): a in S]; 
  AA :=[Eltseq(V*Transpose(Matrix([D])))[1] : V in A]; 
   
//assert &and([Valuation(phi(AA[i])/AA[i]-set[i]) ge Precision(set[1])-Degree(L): i in [1..#set]] );

 //"d is the order of phi which is Inertia degree pf OL and phi =phi'^{d'}";
      //  for j in [3..d] do
       //     Append(~A,c* A[j-1]@phi);
      //  end for;
      //  c_frob := &+(A);
      //  if Norm(c_frob) eq 0 or Valuation(c_frob) ge 1 then
      //     c_frob := Hilbert_Ninety1(c,phi,basis);
     //    end if;
   //     Append(~seq, c_frob);
   // end for;
return AA;
end intrinsic;



intrinsic fix_one(x,p,v,I, iCp)->.{}
    v2 := Valuation(I,p);
    J := I div p^v2;
    m := Order(p);//MaximalOrder(CoefficientField(A));
    x := iCp(x);
    if v lt 0 then
      x *:= Minimum(p)^-v;
    end if;
    v1 := Valuation(x, p);
    y := ChineseRemainderTheorem(J, p^(v1+v2), m!1, m!x);
   /* if #inf ne 0 then //"only when infinite place is there in modulus";
      y := ChineseRemainderTheorem(J*p^(v1+v2), inf, y, [1: t in inf]);
    end if;*/
    assert (x-y) in p^(v1+v2);
    assert (y-1) in J;
    assert Valuation(x/y-1, p) ge v2;
    assert Valuation(y*p^-v1, p) eq 0;
    return y*p^-v1;
  end intrinsic;



intrinsic fix_val(x,p,v,m)->.{}
    if Valuation(x,p) eq 0 then 
       return x;
    end if;
   // v := Minimum([Valuation((x@h), p) : x in gg]);//fix before string
    v2 := Valuation(m,p);
    J := m div p^v2;
    O := Order(p);//MaximalOrder(CoefficientField(A));
    //x := iCp(x);// x in O;
    if v lt 0 then
      x *:= Minimum(p)^-v;
    end if;
    v1 := Valuation(x, p);
    y := ChineseRemainderTheorem(J, p^(v1+v2), O!1, O!x);
   /* if #inf ne 0 then //"only when infinite place is there in modulus";
      y := ChineseRemainderTheorem(J*p^(v1+v2), inf, y, [1: t in inf]);
    end if;*/
    assert (x-y) in p^(v1+v2);
    assert (y-1) in J;
    assert Valuation(x/y-1, p) ge v2;
    assert Valuation(y*p^-v1, p) eq 0;
    return y*p^-v1;
  end intrinsic;







function convert(elt, Mk, M, mo)
    X := Domain(M);
    Z := Domain(Mk);
    phi := Mk(elt);
    aut := InducedAutomorphism(M, phi, mo);
    return Matrix([Eltseq(aut(X.i)) : i in [1..Ngens(X)]]);
    // use InducedAut here!!! XXX and do it in C
end function;
cm := recformat<CohomologyModule, AutomorphismGroup>;


intrinsic CohomologyModule_check(F :: FldAb,psi:Sub := false) -> ModCoho, Map, Map, Map
{The defining ideal group of F as a cohomology module}

  if assigned F`Record and assigned F`Record`CohomologyModule then
    return Explode(F`Record`CohomologyModule);
  end if;

  k := BaseField(F);
/*  g, _, p := AutomorphismGroup(k);
  if Sub cmpne false then
    g := Sub;
  end if;
*/

  g:= Domain(psi);
//  OL := MaximalOrder(Domain(psi(g.1)));
//  p:= hom<g-> Aut(OL)| x:-> hom<OL-> OL |y:-> y@psi(x^-1)  >    >;
 // psi := hom<g-> Codomain(psi)| x:-> psi(x^-1)>;
  p := psi;
  A, mo := NormGroup(F);
  mo := AbsoluteNorm(mo);
  AA := InvariantRepresentation(Domain(A));
  mAA := Coercion(AA, Domain(A));

  inv := AbelianInvariants(AA);
  mats := [ convert(g.i, p, mAA*A, mo) : i in [1..Ngens(g)]];

  C := CohomologyModule(g, inv, mats);
  Zm := RSpace(Integers(), Ngens(AA));
  mp := map<Zm -> AA | x :-> AA!Eltseq(x), y:-> Zm!Eltseq(y)>;
  // p maps the automorphisms group of k onto automorphisms
  // mAA*AA maps between the ideal group of F and the same group
  //      in Smith form (as used in Cohomology module)
  // mp maps between the RSpace from C to the ideal group.     

  if assigned F`Record and Sub cmpeq false then
    F`Record := rec<cm|CohomologyModule := <C, p, mAA*A, mp>,
                       AutomorphismGroup := F`Record`AutomorphismGroup>;
  else
    F`Record := rec<cm|CohomologyModule := <C, p, mAA*A, mp>>;
  end if;
  return C, p, mAA*A, mp;
end intrinsic;




intrinsic AliLocalBrauerGroup_prec(L::FldNum, P::RngOrdIdl, times : autMap := 0, lfc := false) -> Rec
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
    locBr := ali_local_brauer_group_prec(L,P,times, autMap, lfc);
Append(~L`localBrauerGroups, locBr);
    return locBr;
end intrinsic;

intrinsic ali_local_brauer_group_prec(L, P,times, psi, computeLFC)->.{}
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
     m := times*m;// "enlarging the precision";
    vprint Brauer, 1: "lattice precision:" , m;
    vprint Brauer, 2: "theta: ", theta;

    // Localization with enough precision
    vprint Brauer, 1: "compute completion, prec=", 2*m+2;
    vtime Brauer, 1: LP, iota, psiL := completion_with_precision(L, P, psi, m+10);
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






intrinsic gfcUndecomposed_prec(L::FldNum, p0, times ::RngIntElt : psiL := 0) -> ModTupRng, ModCoho, Map, ModTupRngElt, Rec
{ Computes the global fundamental class for a number field L
  in which the prime p0 is undecomposed.
  Optionally one can pass the Galois action on L as map G->Aut(L/Q).
}
    require IsTotallyReal(L) :
            "Just implemented for totally real fields!";
    require #Decomposition(L,p0) eq 1 :
            "Prime must be undecomposed in L!";
    
    t := Cputime();
    
    if psiL cmpeq 0 then
        G,_, psiL := AutomorphismGroup(L);
        psiL := map< Domain(psiL) -> Codomain(psiL) | x:-> psiL(x^(-1)) >;
    else
        G := Domain(psiL);
    end if;
    
    vprint GFC, 1: "compute primes";
    IndentPush();
   primes := {p0} join {x[1] : x in Factorization(Discriminant(RingOfIntegers(L)))};
    /*
    primes:=[x[1] : x in Factorization(Discriminant(RingOfIntegers(L)))];
    choose:=[p: p in primes| #Decomposition(L,p) eq 1][1];
    primes:={p0} join {choose};
    
    */
    subL := Subfields(L);
    for sub in subL do
        F := sub[1];
        CG, m := ClassGroup(F: Bound := BachBound(F));
        S := &cat([  [Ideal(x[1]) : x in Decomposition(F,p)]  : p in primes]);
        
        CGmod, qCG := quo<CG | [s @@ m : s in S]>;
        while #CGmod gt 1 do
            q := Generator(CGmod.1 @@ qCG @m meet Integers());
            S cat:= [Ideal(x[1]) : x in Decomposition(F,q)];
            CGmod, qCG := quo<CG | [s @@ m : s in S]>;
        end while;
        primes := Setseq({Generator(s meet Integers()) : s in S});
    end for;
    //primes := Setseq({Generator(s meet Integers()) : s in S});
    S := &cat([  [Ideal(x[1]) : x in Decomposition(L,p)]  : p in primes]);
    vprint GFC, 1: primes;
    IndentPop();
    
    OL := RingOfIntegers(L);
    // lat := latticeGlob(L : psi := psiL);
    
    
    vprint GFC, 1: "compute S-units and its G-action";
    IndentPush();
    // compute S-Units and G-action
    vtime GFC, 2: US, mUS := SUnitGroup(S);
    GG := [sig : sig in G];
    vtime GFC, 2: sigUS := SUnitAction(mUS, [psiL(sig) : sig in G],S);
    psiUS := map< G -> Aut(US) | sig :-> sigUS[Index(GG,sig)] >;
    IndentPop();
    
    vprint GFC, 1: "Time for set S:", Cputime(t);
    t := Cputime();
    
    vprint GFC, 1: "construct JL";
    IndentPush();
    lst := [];
    LST :=[];
    thetaAll := [];
    // construct JL
    for p in primes do
        vprint GFC, 1: "prime:", p;
        IndentPush();
        
        PL := Factorization(p*OL)[1,1];
        piL := UniformizingElement(PL);
        
        //m := 0;
        //repeat
        //    m := m+1;
        //until &and([b in lat :  b in Generators(PL^m) ]);
        // create lattice for p
        vprint GFC, 2: "compute lattice";
        t := Cputime();
        if RamificationIndex(PL) eq 1 then
            theta := OL!1;
            m := 0;
        else
            theta, m := lattice(PL, piL, psiL);
        end if;
        Append(~thetaAll, OL!theta);
        vprint GFC, 2: "Time:", Cputime(t);
        m := times*m;
        vprint GFC, 2: "compute completion";
	if m lt 50 then
           LP, iotaL := Completion(L, PL : Precision := Max(200,m));
	else
	   LP, iotaL := Completion(L, PL : Precision := Max(200,m));
	end if;
        GP := [g : g in G | &and([  psiL(g)(x) in PL   : x in Generators(PL)]) ];
        GP := sub< G | GP>;
        psiLP := map< GP -> Aut(LP) | g :-> iotaL^(-1) * psiL(g) * iotaL >;
        vprint GFC, 2: "compute module";
       LIST:=[* LP,iotaL,psiLP*];
   ML, psiML, qML := compute_LPmul_modX1(L, PL, psiL, LIST, theta, m);// "increased by 5 times";
  //      ML, psiML, qML := compute_LPmul_modX(L, PL, UniformizingElement(PL), psiL, iotaL, LP, psiLP, theta, m);
        // induce module
        vprint GFC, 2: "compute induced module";
        indML, psiIndML, RL, kappaML, projML := inducedModule(ML, psiML, G);
        diagL := map< US -> indML | x :-> 
            &+([ x @ mUS @ psiL(RL[i]^(-1)) @ iotaL @ qML @ kappaML[i] : i in [1..#kappaML] ]) >;
       projML_L:= hom<indML->L | x:-> (&+[x@projML[i]@@qML :i in [1..#projML]])@@iotaL   >; 
        vprint GFC, 2: "compute cocycle";
        if p ne p0 then
            // trivial cocycle for this
            c2 := map< car<G, G> -> indML | x :-> Zero(indML) >;
        else
            
            // cocycle for p 
            vtime GFC, 2: brGrp := AliLocalBrauerGroup_prec(L, PL, times : autMap := psiL, lfc);
            
            c := TwoCocycle(brGrp`C, brGrp`lfc); 
            C := Group(brGrp`C);
            // c := map< car<C,C> -> brGrp`M | x :-> c(x) @@ brGrp`f1 >;
            // c2 := map< Domain(c) -> Codomain(c) | x :-> c(x[2]^(-1), x[1]^(-1)) >;
            //testCocycle(c2, brGrp`actM );
            
            c2 := map< car<C,C> -> indML | x :-> c(x) @@ brGrp`f1 @@ brGrp`qM @ iotaL @ qML @ kappaML[1] >;
            // pre-compute images
            ll := [x : x in Domain(c2)];
            vtime GFC, 2: llImg := [c2(x) : x in ll];
            c2 := map< Domain(c2) -> Codomain(c2) | x :-> llImg[Index(ll, x)] >;
            
        end if;
       // Append(~lst, [* indML, psiIndML, diagL, c2 *]);
       Append(~lst, [* indML, psiIndML, diagL, c2, projML_L, qML, iotaL,kappaML  *]);
        Append(~LST, [* projML_L, qML,iotaL *]);
       IndentPop();
    end for;
    
    
    // infinite places
    vprint GFC, 1: "modules for infinite places";
    assert &and([ IsReal(inf) : inf in InfinitePlaces(L) ]);
    psiM := map< sub<G | Id(G)> -> Aut(US) | sig :-> hom< US -> US | [US.i : i in [1..#Generators(US)]]> >;
    indML, psiIndML, RL, kappaML, projML := inducedModule(US, psiM, G);
    diagL := map< US -> indML | x :-> 
            &+([ x @ psiUS(RL[i]^(-1)) @ kappaML[i] : i in [1..#kappaML] ]) >;
    c2 := map< car<G, G> -> indML | x :-> Zero(indML) >;
projML_L:= hom<indML->L | x:-> (&+[x@projML[i] :i in [1..#projML]])@mUS   >;

    Append(~lst, [* indML, psiIndML, diagL,  c2, projML_L *] );
    Append(~LST, [* projML_L  *]);
    IndentPop();
    vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();
    
    vprint GFC, 1: "compute idele group of L";
    IndentPush();
    JL, inclJL, projJL := DirectSum([o[1] : o in lst]);
    // precompute projections
    vtime GFC, 2: projJL2 := [ hom< Domain(p) -> Codomain(p) | [ p(JL.i) : i in [1..#Generators(JL)]] >  : p in projJL ];
    // recomputation of projections using injections much faster
    //projJL := [ hom< JL -> lst[k,1] | 
    //    [ Index(seq,i) eq 0 select lst[k,1]!0 else lst[k,1].(Index(seq,i))  : i in [1..#Generators(JL)]] 
    //    >
    //    where seq := [Index(Eltseq(inclJL[k](lst[k,1].i)),1) : i in [1..#Generators(lst[k,1])]]
    //    : k in [1..#lst]];
    
    vtime GFC, 2: actJL := [ hom< JL -> JL |
        [&+( [ JL.j @ projJL[k] @ lst[k,2](sig) @ inclJL[k] : k in [1..#lst]]) : j in [1..#Generators(JL)]]
        > :  sig in GG];
    psiJL := map< G -> Aut(JL) | sig :-> actJL[ Index(GG, sig) ] >;
    
    gamma := map< car<G, G> -> JL | x :-> &+([ x @ lst[i,4] @ inclJL[i]  : i in [1..#lst] ]) >;
    //gfcId := IdentifyTwoCocycle(CohJL, func< x | gamma(x[1],x[2]) @ f1JL >);
    IndentPop();
    vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();
    
    vprint GFC, 1: "compute idele class group of L";
    IndentPush();
    // diagonal embedding of S-units
    embJL := map< US -> JL | x :-> &+([ x @ lst[i,3] @ inclJL[i] : i in [1..#lst]] ) >;
    // factor out S-Units diagonally
    vtime GFC, 2: B := [g @ embJL : g in Generators(US)];
    CL, qCL := quo<JL | B>;
    // time homBasis := [ [CL.i @@ qCL @ psiJL(sig) @ qCL : i in [1..#Generators(CL)]] : sig in GG];
    // psiCL := map< G -> Aut(CL) | sig :-> 
    //     hom< CL -> CL | homBasis[Index(GG, sig)] >
    // >;
    psiCL := map< G -> Aut(CL) | sig :-> Inverse(qCL)*psiJL(sig)*qCL >;
    IndentPop();
    
    vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();
    
    vprint GFC, 1: "compute cohomology of L";
    IndentPush();
    // compute cohomology
    // make right actions
    psiCLr := map< G -> Aut(CL) | g :-> psiCL(g^(-1)) >;
    vtime GFC, 2: CohL := CohomologyModule(G, CL, psiCLr);
    // second cohom. group
    vtime GFC, 2: H2L := CohomologyGroup(CohL,2);
    f1CL := map< CL -> RSpace(Integers(),Dimension(CohL)) |
        x:-> Eltseq(x),
        y:-> CL!Eltseq(y)
    >;
    IndentPop();
    vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();
    
    vprint GFC, 1: "Identify fundamental class:", Cputime(t);
    gammaC := map< car<G, G> -> CL | x :-> x @ gamma @ qCL >;
    gfcId := IdentifyTwoCocycle(CohL, func< x | gammaC(x[1],x[2]) @ f1CL >);
    
    inclUSJL := map< US -> JL | x :-> x @ diagL @ inclJL[#inclJL] >;
//    Req := [* inclJL,projJL, kappaML, projML, lst  *] ;
   // h1:= hom<car<G, G>->L| x:-> x@TwoCoycle(CohL,gfcId)@@f1CL@@comp`qCL@projJL[1]@lst[1,6,1]@@lst[1,5]@@lst[1,7]@@mr>;
   comp := rec<GFCcomp |
        CohL := CohL, f1CL := f1CL, gfcId := gfcId,
        CL := CL, psiCL := psiCL, qCL := qCL,
        primes := primes, US := US, mUS := mUS,
        kappaInf := kappaML, RepInf := RL, inclJL := inclJL,
        inclUSJL := inclUSJL,
        theta := thetaAll >;
   maps:=[];
   CL_L:=[];
   for i in [1..#primes] do
       h:= hom<car<G, G>->L| x:-> x@TwoCocycle(CohL,gfcId)@@f1CL@@comp`qCL@projJL[i]@lst[i,5,1]@@lst[i,6]@@lst[i,7]>;
       h1:= hom<CL->L| x:-> x@@comp`qCL@projJL[i]@lst[i,5,1]@@lst[i,6]@@lst[i,7],
                       y:-> y@lst[i,7] @ lst[i,6] @lst[i,8,1]@ inclJL[i]@comp`qCL    >;
      Append(~maps,h);
      Append(~CL_L,h1);
   end for;

Req := [*LST,lst,projJL,  maps, CL_L, psiJL,inclJL, psiL  *] ;
	
    
return CohL, f1CL, gfcId,comp,Req ;
//return H2L,CohL, f1CL, gfcId, comp;
end intrinsic;
/*

//"from  local side";
//Korrekt so? ^-1 wegen Transformation von psi in "cocycle_lfc" ...
  Mats:=[ Matrix([ Eltseq(mUS(Inverse(mU)(field_map(F,psi(G.j^-1))(mU(Inverse(mUS)(US.i)))))) : i in [1..d]])
        : j in [1..#Generators(G)+1] ];
  //5;
  CM:=CohomologyModule(G, AbelianInvariants(US), Mats);
  H2:=CohomologyGroup(CM,2);//H2;
  //5;
  function Gamma(tup)
    return Eltseq(mUS(Inverse(mU)(gamma(<tup[1],tup[2]>))));
  end function;



*/

intrinsic group_extension_val(L,p ,m)->.
 { group extension for cyclic extensions}
CohL,f1CL,gfc, comp,req:=gfcUndecomposedcl(L,p);
LST := req[1];
lst := req[2];
G:= Group(CohL);
qCL := comp`qCL;
projJL := req[3];
list := [* CohL, gfc, comp *];
B:=[*  *];
u:= TwoCocycle(CohL,gfc);
  for i in [1..#lst] do
        h := hom<car<G,G>->L | x:-> x@u@@f1CL@@qCL@ projJL[i]@lst[i,5]>;
        Append(~B,h);
   end for;
 O := MaximalOrder(L);
 md := m*O;
 r,mr:=RayClassGroup(m*MaximalOrder(L));
 A:= AbelianExtension(mr);
 psi := req[8];
 psi := map< Domain(psi) -> Codomain(psi) | x:-> psi(x^(-1)) >;
 CM,w1,w2,w3 := CohomologyModule_check(A,psi);
 C2:= CohomologyGroup(CM, 2);
mm:= hom<car<G,G>->L| x:-> &*[(<x[1], x[2]>)@h : h in B]>;
// mm:= hom<car<G,G>->L| x:-> &*[x@h : h in B]>;
pp := [x[1]: x in Factorisation( md)];
im_arr := [<X, mm(X)> : X in car<G, G>];
 if  #[x: x in im_arr| x[2] eq 0] gt 0 then 
   error  "remove one projection from CL to L";
 end if;
vv := [Minimum([ Valuation(x[2] ,pp[i]): x in im_arr]): i in [1..#pp]];
vv := [vv[i] div RamificationDegree(pp[i]): i in [1..#pp]];
function fix_one(x)
   aa:=[Valuation(x,p) eq 0 : p in pp[1..#pp]];
   if &and(aa) or x eq 0 then 
      return x;
    else p := [p : p in pp |  Valuation(x,p) ne 0 ][1];
 end if;
  v := vv[Position(pp,p)];
   v2 := Valuation(md, p);
  J := md  div p^v2;
  e := RamificationDegree(p);
  if v lt 0 then
      x *:= Minimum(p)^-v;
    end if;
  v1 := Valuation(x, p);
  y := ChineseRemainderTheorem(J, p^(v1+v2), O!1, O!x);

 assert (x-y) in p^(v1+v2);
    assert (y-1) in J;
    assert Valuation(x/y-1, p) ge v2;
    assert Valuation(y*p^-v1, p) eq 0;
return y*p^-v1;
  end function;
//if  #{x:x in vv} eq 1 and vv[1] eq 0 then
 //  im_arr := [<x[1], x[2]@@mr@@w3> : x in im_arr];
//else  
 im_arr := [<x[1], fix_one(x[2])@@mr@@w3> : x in im_arr];
//end if;
ff := func< X | im_arr[Position([x[1] : x in im_arr], X)][2]>;
 cocycle := IdentifyTwoCocycle(CM, ff);
 //V:=[[<x, x@h>: x in car<G,G>]: h in B];
 //cocycle := IdentifyTwoCocycle(CM,func<x | mm(x[1],x[2])@@mr@@w3>);
 if Type(cocycle) eq ModTupRngElt then
   return Extension(CM, cocycle), CM, cocycle, list;
else return " cocyle is not found, try with different modulus ";
end if;
end intrinsic;



 intrinsic group_extension_prec(L,p,times ,m)->.
 { group extension for cyclic extensions}
CohL,f1CL,gfc, comp,req:=gfcUndecomposed_prec(L,p,times);
LST := req[1];
lst := req[2];
G:= Group(CohL);
qCL := comp`qCL;
projJL := req[3];
list := [* CohL, gfc, comp *];
B:=[*  *];
u:= TwoCocycle(CohL,gfc);
  for i in [1..#lst] do
	h := hom<car<G,G>->L | x:-> x@u@@f1CL@@qCL@ projJL[i]@lst[i,5]>;
	Append(~B,h);
   end for;
 r,mr:=RayClassGroup(m*MaximalOrder(L));   
 A:= AbelianExtension(mr);
 psi := req[8];
 CM,w1,w2,w3 := CohomologyModule_check(A,psi);
 C2:= CohomologyGroup(CM, 2); 
mm:= hom<car<G,G>->L| x:-> &*[(<x[1], x[2]>)@h : h in B]>;
// mm:= hom<car<G,G>->L| x:-> &*[x@h : h in B]>;
 cocycle := IdentifyTwoCocycle(CM,func<x | mm(x[1],x[2])@@mr@@w3>);
 if Type(cocycle) eq ModTupRngElt then 
   return Extension(CM, cocycle), CM, cocycle, list;
else return " cocyle is not found, try with different modulus ";
end if;
end intrinsic;


//////////////////////////////////////////non-cyclic case////////////////////////////




intrinsic gfcCompositum_prec(L::FldNum, L1::FldNum, times) -> ModCoho, Map, ModTupRngElt, Rec
{ Given an arbitrary Galois extension L/Q and a cyclic extension L1/Q
  of the same degree, this method computes then global fundamental
  class of L/Q.
}

   // require IsCyclic(L1) and Degree(L) eq Degree(L1) :
     //       "Second number field must be cyclic and of the same degree!";
require IsCyclic(L1):"Second number field must be cyclic";
    t := Cputime();

    vprint GFC, 1: "compute composite field";
    IndentPush();
    vtime GFC, 1: N := OptimizedRepresentation(Compositum(L,L1));
    assert IsTotallyReal(N);
    ON := RingOfIntegers(N);

    Gamma, _ ,psiN := AutomorphismGroup(N);
    psiN := map< Domain(psiN) -> Codomain(psiN) | x :-> psiN(x^(-1)) >;
    IndentPop();

    OL := RingOfIntegers(L);

    vprint GFC, 1: "compute primes";
    IndentPush();
    
    primes:=[];
primes := [f[1] : f in Factorization(Discriminant(ON))];
/*  seq := [p :   p in primes | #Decomposition(L1, p) eq 1];
    if #seq eq 0 then
    p0:=findUndecomposedPrime(L1);
    primes:= Sort([p0] cat primes);
    else p0:=Sort(seq)[1];
    end if;*/
    vtime GFC, 2: primes := trivialSClassNumberPrimes_check(L,L1,N : primes := primes);
   // prime:=trivialSclassless(L,L1,N);
   // primes:=&cat[prime,primes];
   // set:={x: x in primes};
   // primes:=[x : x in set];
     seq := [p :   p in primes | #Decomposition(L1, p) eq 1];
if #seq eq 0 then
    p0:=findUndecomposedPrime(L1);
    primes:= Sort([p0] cat primes);
    else p0:=Sort(seq)[1];
    end if;
if #seq gt 1 and seq[2] lt 50 then
   p0 := Sort(seq)[2];;
end if;



// vtime GFC, 2: primes := trivialSClassNumberPrimes(N : primes := primes);
    S := &cat([  [Ideal(x[1]) : x in Decomposition(N,p)]  : p in primes]);
    vprint GFC, 1: primes;
    IndentPop();

    vprint GFC, 1: "compute S-units and its G-action";
    IndentPush();
    // compute S-Units and G-action
    vtime GFC, 1: US, mUS := SUnitGroup(S);
    GammaSeq := [sig : sig in Gamma];
    vtime GFC, 1: sigUS := SUnitAction(mUS, [psiN(sig) : sig in GammaSeq],S);
    psiUS := map< Gamma -> Aut(US) | sig :-> sigUS[Index(GammaSeq,sig)] >;
    // S-units for L
    //H := FixedGroup(N,L);
    //K:=[ Kernel(VerticalJoin(Matrix([  Eltseq(US.i @ SUnitAction(mUS, psiN(h),S) - US.i)  :  i in [1..#Generators(US)]]),Transpose(D))) : h in H ] where D := DiagonalMatrix([Order(US.i) :  i in [1..#Generators(US)] ]);
    //K := [ Kernel(Transpose(HorizontalJoin(
    //Transpose(Matrix([  Eltseq(US.i @ SUnitAction(mUS, psiN(h),S) - US.i)  :  i in [1..#Generators(US)]])), D)))
    //: h in H ]
    //where D := DiagonalMatrix([Order(US.i) :  i in [1..#Generators(US)] ]); "is same as the next";
   
   
   
   
   H := FixedGroup(N,L);
   H1 := FixedGroup(N,L1);
    K := [ Kernel(Transpose(HorizontalJoin(
        Transpose(Matrix([  Eltseq(US.i @ psiUS(h) - US.i)  :  i in [1..#Generators(US)]])), D)))
        : h in H ]
        where D := DiagonalMatrix([Order(US.i) :  i in [1..#Generators(US)] ]);
    USL := &meet([sub<US| [US!Eltseq(b)[1..#Generators(US)] :  b in Basis(k)]> : k in K]);

 assert &and([ g @ mUS in L : g in Generators(USL)]);
    IndentPop();

    vprint GFC, 1: "Time for set S:", Cputime(t);
    t := Cputime();

    vprint GFC, 1: "construct JN";
    IndentPush();
    lst := [];
    LST := [];
    thetaAll := [];
    maps := [];
    for p in primes do
        vprint GFC, 1: "prime:", p;
        IndentPush();
        PN := Factorization(p*ON)[1,1];
        piN := UniformizingElement(PN);

        vprint GFC, 2: "compute lattice";
        t := Cputime();
        if RamificationIndex(PN) eq 1 then
            theta := ON!1;
            m := 0;
        else
            theta, m := lattice(PN, piN, psiN);
            for i in [1..2] do
                theta1, m1 := lattice(PN, piN, psiN);
                if m1 lt m then
                    theta := theta1;
                    m := m1;
                end if;
            end for;
        end if;
        Append(~thetaAll, ON!theta);
        vprint GFC, 2: "Time:", Cputime(t);

        /*
        print "compute completion, prec =", Max(100,m*2);
        NP, iotaN := Completion(N, PN : Precision := Max(100,m*2));
        GammaP := [g : g in Gamma | &and([  psiN(g)(x) in PN   : x in Generators(PN)]) ];
        GammaP := sub< Gamma | GammaP>;
        psiNP := map< GammaP -> Aut(NP) | g :-> iotaN^(-1) * psiN(g) * iotaN >;
        */
        //print "completion with sufficient precicion for computations up to precision ", m+10;
       m := times*m;
	 vprint GFC, 2: "compute completion, prec =", m+10;
       if p eq 2 then 
         NP, iotaN := Completion(N, PN : Precision := Max(200,m*2));
        GammaP := [g : g in Gamma | &and([  psiN(g)(x) in PN   : x in Generators(PN)]) ];
        GammaP := sub< Gamma | GammaP>;
        psiNP := map< GammaP -> Aut(NP) | g :-> iotaN^(-1) * psiN(g) * iotaN >;

	else NP, iotaN, psiNP := completion_with_precision(N,PN,psiN, Max(200,m+10));
	 
      // vtime GFC, 2: NP, iotaN, psiNP := completion_with_precision(N,PN,psiN, m+10);
       end if;
        LIST:=[*NP,iotaN,psiNP*];
        GammaP := Domain(psiNP);
        vprint GFC, 2: "compute module";
       // if p eq 2 then
//	   MN, psiMN, qMN := compute_LPmul_modX_check(N, PN, psiN,LIST, theta, m);
	//else
        vtime GFC, 2: MN, psiMN, qMN := compute_LPmul_modX1(N, PN, psiN,LIST, theta, m);
	//end if;
//    vtime GFC, 2: MN, psiMN, qMN := compute_LPmul_modX(N, PN, piN, psiN, iotaN, NP, psiNP, theta, m);  
	// induce module
        vprint GFC, 2: "compute induced module";
       H := FixedGroup(N,L);
        //R := [Gamma!x : x in r] where r := leftCosetRepresentatives(H, H meet GammaP);
        indMN, psiIndMN, RN, kappaMN, projMN := inducedModule(MN, psiMN, Gamma);// : RepSys := R);
        diagN := map< N -> indMN | x :->
            &+([ x @ psiN(RN[i]^(-1)) @ iotaN @ qMN @ kappaMN[i] : i in [1..#kappaMN] ]) >;

// H := FixedGroup(N,L);
        K := [ Kernel(Transpose(HorizontalJoin(
            Transpose(Matrix([  Eltseq(indMN.i @ psiIndMN(h) - indMN.i)  :  i in [1..#Generators(indMN)]])), D)))
        : h in H ]
        where D := DiagonalMatrix([Order(indMN.i) :  i in [1..#Generators(indMN)] ]);
        indML := &meet([sub<indMN | [indMN!Eltseq(b)[1..#Generators(indMN)] :  b in Basis(k)]> : k in K]);

        assert (N!L.1) @ diagN in indML;
        /*
        H := FixedGroup(N,L1);
        K := [ Kernel(Transpose(HorizontalJoin(
            Transpose(Matrix([  Eltseq(indMN.i @ psiIndMN(h) - indMN.i)  :  i in [1..#Generators(indMN)]])), D)))
            : h in H ]
            where D := DiagonalMatrix([Order(indMN.i) :  i in [1..#Generators(indMN)] ]);
        indML1 := &meet([sub<indMN | [indMN!Eltseq(b)[1..#Generators(indMN)] :  b in Basis(k)]> : k in K]);
        assert (N!L1.1) @ diagN in indML1;
        */

        if p ne p0 then
            // trivial cocycle for this
            c2 := map< car<Gamma, Gamma> -> indMN | x :-> Zero(indMN) >;
        else
            vprint GFC, 2: "compute cocycle, prec =", m;
            // compute cocycle for p
        //    H := FixedGroup(N,L1);
            C, qC := quo< Gamma | H1>;
            //psiL1 := map< C -> Aut(L1) | g :-> Coercion(L1,N) * psiN(g @@ qC) * Coercion(N,L1) >;
            psiL1 := map< C -> Aut(L1) | g :->
                hom< L1 -> L1 | L1.1 @ Coercion(L1,N) @ psiN(g @@ qC) @ Coercion(N,L1) >
            >;

            // compute ML1
            K := [ Kernel(Transpose(HorizontalJoin(
                Transpose(Matrix([  Eltseq(indMN.i @ psiIndMN(h) - indMN.i)  :  i in [1..#Generators(indMN)]])), D)))
                : h in H1 ]
                where D := DiagonalMatrix([Order(indMN.i) :  i in [1..#Generators(indMN)] ]);
            indML1 := &meet([sub<indMN | [indMN!Eltseq(b)[1..#Generators(indMN)] :  b in Basis(k)]> : k in K]);
            psiIndML1 := map< C -> Aut(indML1) |
                sig :-> Coercion(indML1, indMN)*psiIndMN(sig @@ qC)*Coercion(indMN,indML1) >;

            // compute completion of L1
            PL1 := Factorization(p*RingOfIntegers(L1))[1,1];
            //print "completion with sufficient precicion for computations up to precision ", m+10;
            vprint GFC, 2: "compute completion, prec =", m+10;
            L1P, iotaL1, psiL1P := completion_with_prec(L1,PL1,psiL1, Max(200,m+20));
           // L1P, iotaL1 := Completion(L1, PL1 : Precision := Max(200,5*m )); //Max(100,m*2));
            //psiL1P := map< C -> Aut(L1P) | g :-> iotaL1^(-1) * psiL1(g) * iotaL1 >;
            // cocycle C x C -> L1P
            //SetVerbose("CocycleLFC", 1);
            //  c := ClCocycleLFC(L1P, pAdicField(L1P), m+5 : psi := psiL1P);
	    if p gt 50 then
             c := CLocalFundamentalClassSerre_check(L1P, pAdicField(L1P), m+5 : psi := psiL1P);
            else c := ClCocycleLFC(L1P, pAdicField(L1P), m+5 : psi := psiL1P);
	    //c :=CLocalFundamentalClassSerre(L1P, pAdicField(L1P), m+5 : psi := psiL1P);
	    end if;
	    // inflation
            c2 := map< car<Gamma,Gamma> -> indMN | x :-> c(x[1]@qC, x[2]@qC) @@ iotaL1 @ Coercion(L1,N) @ diagN>;
            vprint GFC, 2: "test cocycle";
            vtime GFC, 2: assert testCocycleGenerators(c2, psiIndMN );
            c2 := map< Domain(c2) -> Codomain(c2) | x:-> c2(x[2]^(-1), x[1]^(-1)) >;
 end if;

        diagN := mUS*diagN;
        Append(~lst, [* indML, indMN, psiIndMN, diagN, c2 *]);
        Append(~LST, [* PN,m,RN,iotaN, qMN,kappaMN, projMN *]);
	Append(~maps, [* projMN, qMN, indML*]);
        IndentPop();
    end for;

  // infinite places
    vprint GFC, 1: "modules for infinite places";
    assert &and([ IsReal(inf) : inf in InfinitePlaces(N) ]);
    psiM := map< sub<Gamma | Id(Gamma)> -> Aut(US) | sig :-> hom< US -> US | [US.i : i in [1..#Generators(US)]]> >;
    indMN, psiIndMN, RN, kappaMN, projMN := inducedModule(US, psiM, Gamma);
    diagN := map< US -> indMN | x :->
            &+([ x @ psiUS(RN[i]^(-1)) @ kappaMN[i] : i in [1..#kappaMN] ]) >;
    c2 := map< car<Gamma, Gamma> -> indMN | x :-> Zero(indMN) >;
    // Fix-module by H
    H := FixedGroup(N,L);
    K := [ Kernel(Transpose(HorizontalJoin(
        Transpose(Matrix([  Eltseq(indMN.i @ psiIndMN(h) - indMN.i)  :  i in [1..#Generators(indMN)]])), D)))
        : h in H ]
        where D := DiagonalMatrix([Order(indMN.i) :  i in [1..#Generators(indMN)] ]);
    indML := &meet([sub<indMN | [indMN!Eltseq(b)[1..#Generators(indMN)] :  b in Basis(k)]> : k in K]);
    assert &and([ x @ diagN in indML : x in Generators(USL)]);

    Append(~lst, [* indML, indMN, psiIndMN, diagN, c2 *] );
    Append(~LST, [* PN, RN,psiIndMN, kappaMN, projMN, psiUS *]);
    IndentPop();

vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();

    // Finitely generated idele group
    vprint GFC, 1: "compute idele group of N";
    IndentPush();
    JN, inclJN, projJN := DirectSum([o[2] : o in lst]);
    // recompute projections
    vtime GFC, 1: projJN := [ hom< Domain(p) -> Codomain(p) | [ p(JN.i) : i in [1..#Generators(JN)]] >  : p in projJN ];
    //projJN := [ hom< JN -> lst[k,2] |
       // [ Index(seq,i) eq 0 select lst[k,2]!0 else lst[k,2].(Index(seq,i))  : i in [1..#Generators(JN)]]
      //  >
    //    where seq := [Index(Eltseq(inclJN[k](lst[k,2].i)),1) : i in [1..#Generators(lst[k,2])]]
  //      : k in [1..#lst]];

    vtime GFC, 1: actJN := [ hom< JN -> JN |
        [&+( [ JN.j @ projJN[k] @ lst[k,3](sig) @ inclJN[k] : k in [1..#lst]]) : j in [1..#Generators(JN)]]
        > :  sig in GammaSeq];
    psiJN := map< Gamma -> Aut(JN) | sig :-> actJN[ Index(GammaSeq, sig) ] >;

    gamma := map< car<Gamma, Gamma> -> JN | x :-> &+([ x @ lst[i,5] @ inclJN[i]  : i in [1..#lst] ]) >;
    //gammaL := map< Domain(gamma) -> Codomain(gamma) | x :-> gamma(x[2]^(-1), x[1]^(-1)) >;
    //time testCocycleGenerators(gammaL, psiJN);
    IndentPop();

 vprint GFC, 1: "Time:", Cputime(t);
    t := Cputime();

    vprint GFC, 1: "compute idele class group of N";
    IndentPush();
    // diagonal embedding of S-units
    embJN := map< US -> JN | x :-> &+([ x @ lst[i,4] @ inclJN[i] : i in [1..#lst]] ) >;
    // factor out S-Units diagonally
    vtime GFC, 1: B := [g @ embJN : g in Generators(US)];
    CN, qCN := quo<JN | B>;
    psiCN := map< Gamma -> Aut(CN) | sig :-> Inverse(qCN)*psiJN(sig)*qCN >;
    //gammaL := map< Domain(gamma) -> CN | x :-> gamma(x[2]^(-1), x[1]^(-1)) @ qCN >;
    //time testCocycleGenerators(gammaL, psiCN);
    IndentPop();
//"till here works well in 260 seconds around for S3":

 vprint GFC, 1: "compute cohomology of N";
    IndentPush();
    // compute cohomology
    // make right actions
    psiCNr := map< Gamma -> Aut(CN) | g :-> psiCN(g^(-1)) >;
    vtime GFC, 1: CohN := CohomologyModule(Gamma, CN, psiCNr);
    f1CN := map< CN -> RSpace(Integers(),Dimension(CohN)) | x:-> Eltseq(x), y:-> CN!Eltseq(y) >;
    // second cohom. group
    //time H2N := CohomologyGroup(CohN,2);
    vtime GFC, 1: H1N := CohomologyGroup(CohN,1);

    gammaC := map< car<Gamma, Gamma> -> CN | x :-> x @ gamma @ qCN >;
    //gfcId := IdentifyTwoCocycle(CohN, func< x | gammaC(x[1],x[2]) @ f1CN >);
    IndentPop();
//"till here also works well in 360 seconds around for S3":    

vprint GFC, 1: "Time for cohomology of N:", Cputime(t);
    t := Cputime();

    vprint GFC, 1: "compute idele group of L";
    // Cohomology of L
    JL, inclJL, projJL := DirectSum([o[1] : o in lst]);

    embLN := map< JL -> JN |
        x :-> &+([ x @ projJL[i] @ Coercion(lst[i,1], lst[i,2]) @ inclJN[i] : i in [1..#lst]]),
        y :-> &+([ y @ projJN[i] @ Coercion(lst[i,2], lst[i,1]) @ inclJL[i] : i in [1..#lst]])
    >;
    G, qG := quo< Gamma | FixedGroup(N,L) >;
    psiJL := map< G -> Aut(JL) | sig :-> embLN * (sig @@ qG @ psiJN) * Inverse(embLN) >;

     vprint GFC, 1: "compute idele class group of L";
    IndentPush();
    vtime GFC, 1: B := [g @ embJN @@ embLN : g in Generators(USL)];
    CL, qCL := quo<JL | B>;

 psiCL := map< G -> Aut(CL) | sig :-> Inverse(qCL)*psiJL(sig)*qCL >;
    IndentPop();

    vprint GFC, 1: "compute cohomology of L";
    IndentPush();
    // compute cohomology
    // make right actions
    psiCLr := map< G -> Aut(CL) | g :-> psiCL(g^(-1)) >;
    vtime GFC, 1: CohL := CohomologyModule(G, CL, psiCLr);
    // second cohom. group
    vtime GFC, 1: H2L := CohomologyGroup(CohL,2);
   assert #H2L eq Degree(L);
    f1CL := map< CL -> RSpace(Integers(),Dimension(CohL)) | x:-> Eltseq(x), y:-> CL!Eltseq(y) >;
    IndentPop();


    vprint GFC, 1: "Time for all the computation for L:", Cputime(t);
    t := Cputime();
     psiL := map< G -> Aut(L) | g :-> Coercion(L,N) * psiN(g @@ qG) * Coercion(N,L) >;
    mUSL := map< USL -> L | x :-> L!(x @ mUS) >;
    inclUSJL := map< USL -> JL | x :-> (US!x) @ diagN @ inclJL[#inclJL] >;

    comp := rec<GFCcomp |
        CohL := CohL, f1CL := f1CL, //gfcId := gfcId,
        CL := CL, psiCL := psiCL, qCL := qCL,
        primes := primes, US:= USL, mUS := mUSL,
        //kappaInf := kappaML, RepInf := RL, inclJL := inclJL,
        inclUSJL := inclUSJL,
        theta := thetaAll >;

    Req:=[* LST,lst, projJL, inclJL,qG,gammaC,qCN,psiL,embLN,CohN,f1CN *];
 
vprint GFC, 1: "find global fundamental class of L";
    IndentPush();
   for k in [ i : i in [1..Degree(L)] | GCD(i,Degree(L)) eq 1 ] do
    // for k in [ i : i in [1..Degree(L)]] do
     vprintf GFC, 1: ".";
        c := TwoCocycle(CohL, k*H2L.1);
        c2 := map< car< G, G> -> CL | x :-> c(<x[1],x[2]>) @@ f1CL >;
        c3 := map< car< Gamma, Gamma> -> CN | x :-> c2(x[1]@qG,x[2]@qG) @@ qCL @ embLN @ qCN>;
        //c4 := func< x | c3(x) @ f1CN>;
        dif := map< Domain(c3) -> Codomain(c3) | g :-> gammaC(g)-c3(g) >;
        bool, prog := IsTwoCoboundary(CohN, func< x | dif(x[1],x[2]) @ f1CN >);
        if bool then
            vprint GFC, 1: " found.";
            IndentPop();
            comp`gfcId := k*H2L.1;
           return CohL, f1CL, k*H2L.1, comp,Req;
        end if;
    end for;
    vprint GFC, 1: " failed.";
    IndentPop();
    error "Global fundamental class could not be found!!";
end intrinsic;

 intrinsic group_extension_compositum_val(L,L1, m0)->.
 { construct an extension with fundamental class for any abelian extension}
// *"Construct>>>>>>psiindML>>>> with each prime";
 CohL, f1CL,gfc,comp, req := gfcCompositumcl(L,L1);
 list := [* CohL, f1CL,gfc,comp, req *];
G:=Group(CohL);
 CL:=comp`CL;
qCL:=comp`qCL;
 LST:=req[1];
lst:=req[2];
 u:= TwoCocycle(CohL,gfc);
projJL:= req[3];
//psi := req[8];
mUSL:= comp`mUS;
projMN:= LST[#lst,5]; 

Hom :=[* *];
for i in [1..#LST-1] do
  P := LST[i,1];
  m := LST[i,2];
  map :=  LST[i,4]*LST[i,5];
  Append(~Hom, solve_reconstruction(P,m,L,map));
end for;  
 
//h2:=hom<car<G,G>->L |  x:->reconstruction(LST[i,1],  &+[x@u@@f1CL@@qCL@projJL[i]@ Coercion( lst[i,1], lst[i,2] )@ lst[i,3](LST[i,3][k])@ LST[i,7][k]: k in [1..#LST[i,7]]]@@LST[i,5]@@LST[i,4], LST[i,2],L)> where i :=2;

set := [hom<car<G,G>->L |  x:-> (&+ [x@u@@f1CL@@qCL@projJL[i]@ Coercion( lst[i,1], lst[i,2] )@ lst[i,3](LST[i,3][k])@ LST[i,7][k]: k in [1..#LST[i,7]]]) @Hom[i] >: i in [1..#lst-1]];

h_inf:= hom<car<G,G>->L |  x:-> (&+[(x@u@@f1CL@@qCL@projJL[i]@Coercion(lst[i,1],lst[i,2]))@projMN[j]: j in [1..#projMN]]) @mUSL> where i :=#lst;
//Append(~set,h_inf);
r,mr := RayClassGroup(m0*MaximalOrder(L));
m0 := m0 * MaximalOrder(L);
//B:=[* Hom[1], Hom[2], h_inf *];
A:=AbelianExtension(mr);
//ca,w1,w2,w3 := CohomologyModule_check(A,psi);
ca,w1,w2,w3 := CohomologyModule(A);
CohomologyGroup(ca,2);
H := Group(ca);
_,phi:= IsIsomorphic(H,G);
 mm:= hom<car<G,G>-> L | x:->&*[x@h:h in set ]>;
cocycle :=IdentifyTwoCocycle(ca, func<x | fix_val_mod((<x[1]@phi, x[2]@phi>)@mm, m0)@@mr@@w3>);
if Type(cocycle) eq ModTupRngElt then
   return Extension(ca, cocycle), ca,cocycle,list;
 else return "no cocycle,try with other modulus"  ;
end if ;

end intrinsic;


intrinsic fix_val_mod(x,m)->.{}

  pp:= [x[1]: x in Factorisation(m)];

  vv := [Valuation(x,p): p in pp];
  if &and[a eq 0 : a in vv] then
     return x ;
  else
       num := Position(vv,[a: a in vv | a ne 0][1]);
       p:= pp[num ];
  end if;
   v := Valuation(x,p);

   /* if Valuation(x,p) eq 0 then
       return x;
    end if;*/
   // v := Minimum([Valuation((x@h), p) : x in gg]);//fix before string
    v2 := Valuation(m,p);
    J := m div p^v2;
    O := Order(p);//MaximalOrder(CoefficientField(A));
    //x := iCp(x);// x in O;
    if v lt 0 then
      x *:= Minimum(p)^-v;
    end if;
    v1 := Valuation(x, p);
    y := ChineseRemainderTheorem(J, p^(v1+v2), O!1, O!x);
   /* if #inf ne 0 then //"only when infinite place is there in modulus";
      y := ChineseRemainderTheorem(J*p^(v1+v2), inf, y, [1: t in inf]);
    end if;*/
    assert (x-y) in p^(v1+v2);
    assert (y-1) in J;
    assert Valuation(x/y-1, p) ge v2;
    assert Valuation(y*p^-v1, p) eq 0;
    return y*p^-v1;
  end intrinsic;



intrinsic group_extension_inv(L,p ,m,inv :GFC:=[])->.
 { group extension for cyclic extensions}
 O := MaximalOrder(L);
 md := m*O;
 r,mr:=RayClassGroup(m*MaximalOrder(L));
 s:= Subgroups(r: Quot:= inv);
if #s eq 0 then 
   error "No subgroup of given invariant";
else
   assert #s ge 1;
   AA := [];
   mqq := [];
   for i in [1..#s] do
     sub := s[i];
     _,mq := quo< r| sub`subgroup>;
     Abel := AbelianExtension(Inverse(mq)*mr );
     if IsNormal(Abel : All) then
        Append(~AA, Abel);
        Append(~mqq, mq);
         break i;
     end if;
   end for;
end if;
if #AA eq 0 then
   return "No normal extension of the given data";
end if;
mq := mqq[1];
A := AA[1];
_,m0,_ := NormGroup(A);
 if #GFC ge 1 then
   CohL,f1CL,gfc, comp,req:= Explode(GFC);
 else 
    CohL,f1CL,gfc, comp,req:=gfcUndecomposedcl(L,p);
  end if;
LST := req[1];
lst := req[2];
G:= Group(CohL);
qCL := comp`qCL;
projJL := req[3];
list := [* CohL,f1CL, gfc, comp, req *];
B:=[*  *];
u:= TwoCocycle(CohL,gfc);
  for i in [1..#lst] do
        h := hom<car<G,G>->L | x:-> x@u@@f1CL@@qCL@ projJL[i]@lst[i,5]>;
        Append(~B,h);
   end for;
// A:= AbelianExtension(mr);
 psi := req[8];
 CM,w1,w2,w3 := CohomologyModule_check(A,psi);
//CM,w1,w2,w3 := CohomologyModule(A);
C2:= CohomologyGroup(CM, 2);
mm:= hom<car<G,G>->L| x:-> &*[(<x[1], x[2]>)@h : h in B]>;
pp := [x[1]: x in Factorisation( md)];
im_arr := [<X, mm(X)> : X in car<G, G>];
 if  #[x: x in im_arr| x[2] eq 0] gt 0 then
   error  "remove one projection from CL to L";
 end if;
vv := [Minimum([ Valuation(x[2] ,pp[i]): x in im_arr]): i in [1..#pp]];
vv := [vv[i] div RamificationDegree(pp[i]): i in [1..#pp]];
function fix_one(x)
   aa:=[Valuation(x,p) eq 0 : p in pp[1..#pp]];
   if &and(aa) or x eq 0 then
      return x;
    else p := [p : p in pp |  Valuation(x,p) ne 0 ][1];
 end if;
  v := vv[Position(pp,p)];
   v2 := Valuation(md, p);
  J := md  div p^v2;
  e := RamificationDegree(p);
  if v lt 0 then
      x *:= Minimum(p)^-v;
    end if;
  v1 := Valuation(x, p);
  y := ChineseRemainderTheorem(J, p^(v1+v2), O!1, O!x);

/* if #inf ne 0 then //"for infinite place"
      y := ChineseRemainderTheorem(J*p^(v1+v2), inf, y, [1: t in inf]);
    end if;  

*/
 assert (x-y) in p^(v1+v2);
    assert (y-1) in J;
    assert Valuation(x/y-1, p) ge v2;
    assert Valuation(y*p^-v1, p) eq 0;
return y*p^-v1;
  end function;
 im_arr := [<x[1], fix_one(x[2])@@mr@mq@@w3> : x in im_arr];
//end if;
ff := func< X | im_arr[Position([x[1] : x in im_arr], X)][2]>;
 cocycle := IdentifyTwoCocycle(CM, ff);


//V:=[[<x, x@h>: x in car<G,G>]: h in B];
 //cocycle := IdentifyTwoCocycle(CM,func<x | mm(x[1],x[2])@@mr@@w3>);
 if Type(cocycle) eq ModTupRngElt then
   return Extension(CM, cocycle), CM, cocycle, list;
else return " cocyle is not found, try with different modulus ";
end if;
end intrinsic;




intrinsic group_extension_compositum_inv(L,L1,m0, inv:GFC:=[])->.
 { construct an extension with fundamental class for any abelian extension}
// *"Construct>>>>>>psiindML>>>> with each prime";
if #GFC ge 1 then
 CohL, f1CL,gfc,comp, req := Explode(GFC);
else
 time CohL, f1CL,gfc,comp, req := gfcCompositumcl(L,L1);
end if;
 list := [* CohL, f1CL,gfc,comp, req *];
G:=Group(CohL);
 CL:=comp`CL;
qCL:=comp`qCL;
 LST:=req[1];
lst:=req[2];
 u:= TwoCocycle(CohL,gfc);
projJL:= req[3];

mUSL:= comp`mUS;
projMN:= LST[#lst,5];
primes := comp`primes;
Hom :=[* *];
for i in [1..#LST-1] do
  P := LST[i,1];
  m := LST[i,2];
  map :=  LST[i,4]*LST[i,5];
  Append(~Hom, solve_reconstruction(P,m,L,map));
end for;

set := [hom<car<G,G>->L |  x:-> (&+ [x@u@@f1CL@@qCL@projJL[i]@ Coercion( lst[i,1], lst[i,2] )@ lst[i,3](LST[i,3][k])@ LST[i,7][k]: k in [1..#LST[i,7]]]) @Hom[i] >: i in [1..#lst-1]];

h_inf:= hom<car<G,G>->L |  x:-> (&+[(x@u@@f1CL@@qCL@projJL[i]@Coercion(lst[i,1],lst[i,2]))@projMN[j]: j in [1..#projMN]]) @mUSL> where i :=#lst;
//Append(~set,h_inf);
r,mr := RayClassGroup(m0*MaximalOrder(L));
s:= Subgroups(r: Quot:= inv);

assert #s ge 1;
AA := [];
mqq := [];
for i in [1..#s] do
    sub := s[i];
  _,mq := quo< r| sub`subgroup>;
  Abel := AbelianExtension(Inverse(mq)*mr );
  if IsNormal(Abel : All) then
       Append(~AA, Abel);
       Append(~mqq, mq);
       break i;
    end if;
end for;
if #AA eq 0 then
   return "No normal extension of the given data";
end if;
mq := mqq[1];
A := AA[1];
_,m0,_ := NormGroup(A);

ca,w1,w2,w3 := CohomologyModule(A);
CohomologyGroup(ca,2);
H := Group(ca);
_,phi:= IsIsomorphic(H,G);
// mm:= hom<car<G,G>-> L | x:->&*[x@h:h in set ]>;
  mm_fin:= hom<car<G,G>-> r | x:->&+[ fix_val_mod_S( x@set[i],m0, [primes[i]])@@mr  :  i in [1..#set]]>;
mm_inf := hom< car<G,G>-> r | x:-> fix_val_mod_S(x@h_inf,m0, primes)@@mr  >;
set := [* mm_fin, mm_inf*];
mm := hom< car<G,G> -> r |x:-> &+[x@h : h in set ] >;
cocycle :=IdentifyTwoCocycle(ca, func<x | (<x[1]@phi, x[2]@phi>)  @mm@mq @@w3>);
//cocycle :=IdentifyTwoCocycle(ca, func<x | (<x[1]@phi, x[2]@phi>)@mm@@mr@@w3>);
//cocycle :=IdentifyTwoCocycle(ca, func<x | fix_val_mod_S((<x[1]@phi, x[2]@phi>)@mm, m0, primes)@@mr@@w3>);
if Type(cocycle) eq ModTupRngElt then
   return Extension(ca, cocycle), ca,cocycle,list;
 else return "no cocycle,try with other modulus"  ;
end if ;

end intrinsic;






