

//import "galois.m": eval_field_automorphism, field_map;
declare verbose Segment, 1;
declare verbose group, 1;
// local Brauer group
GaloisGrp := recformat<
   // L    : FldNum,        // the corresponding global field
   // P    : RngOrdIdl,     // with prime ideal
    p    : RngIntElt,
    //m    : RngIntElt      // precision needed

    gamma : Map,        // the module
    E : SeqEnum,          // group action on M
    S : GrpAb,
    Ext_F : GrpFP
   // qM   : Map,          // projection map L^\times ->> M
  //    lfc  : ModTupRngElt  // local fundamental class
>;



///////////////factorisation will be mostly over totally ramified extension so divide the polynomial in smaller degree and then check factorisation over L1/////////////////////////////////////////////
intrinsic eval_field_automorphism(x, m)->.{}
    F := Parent(x);
    R := RingOfIntegers(F);
    pi := UniformizingElement(F);

    xv := Valuation(x);
    xu := ShiftValuation(x, -xv);
    return (F ! m(UniformizingElement(R)))^xv * F ! m(R ! xu);
end intrinsic;

intrinsic field_map(L, m)->.{}
    res := map<L -> L | x :-> eval_field_automorphism(x, m)>;
    res`LocalRingMap := m;
    return res;
end intrinsic;

//////////////////////pReduction for e_0 ne 1////////////////////
function degree_zeta_e(e,q)
  f:=1;
  while not ((q^f-1) mod e) eq 0 do
    f:=f+1;
  end while;
  return f;
end function;


function composite(list,K)
        N:=list[1];
        for i in [2..#list] do
                _,_,E:=Factorization(Polynomial(N,DefiningPolynomial(list[i],K)):Extensions:=true);
                N:=E[1]`Extension;
        end for;
        return N;
end function;


//vgl. Algorithmus 5.1
intrinsic pReduction_e0(f::RngUPolElt) -> FldPad
  {f\in K[x] eisenstein mit p teilt Grad(f) -> Erweiterung T/K, so dass N/T p-Erweiterung ist}

  K:=CoefficientRing(f);
  p:=Prime(K);
  q:=#ResidueClassField(Integers(K));
  n:=Degree(f);
  e0:=n div p^Valuation(n,p);
  flaches_Segment:=false;
  if e0 ne 1 then
    flaches_Segment:=true;
  end if;

 L:=ext<K|f>;
  vpInfos:=XFlatRamificationPolygon(f); // [ ..., < h_i, e_i, p^{s_{i-1}}, A_i(y)>, ... ]
  if flaches_Segment then
    ell:=#vpInfos-1;    //letztes Segment separat behandeln
  else
    ell:=#vpInfos;
  end if;

  b:=[];
  v:=[];
  eff:=[];
  for i in [1..ell] do          //letztes Segment separat behandeln!
    _,bi,_:=XGCD(vpInfos[i][1],vpInfos[i][2]);
    Append(~b,bi);
    vi:=bi*(n div vpInfos[i][3])+n+1;
    Append(~v,vi);
    f1:=LCM([Degree(s[1]) : s in Factorization(vpInfos[i][4])]); //ass. Trägheit
    f0:=degree_zeta_e(e0*vpInfos[i][2],q);
    Append(~eff,LCM(f0,f1));
  end for;

 lcm:=LCM(eff);
  if lcm gt 1 then
    U:=UnramifiedExtension(K,LCM(eff));
  else
    U:=K;
  end if;
  RU,mRU:=ResidueClassField(Integers(U));
  P<x>:=PolynomialRing(U);

  Tis:=[];
  for i in [1..ell] do
    eps:=U!Roots(Polynomial(RU,vpInfos[i][4]))[1][1];
//    Ti:=ext<U|x^(e0*vpInfos[i][2])-((-1)^v[i]*eps^(b[i]*n)*Coefficient(f,0))>;
       Ti:=ext<U|x^(vpInfos[i][2])-((-1)^v[i]*eps^(b[i]*n)*Coefficient(f,0))>;
       // "removing the first segment which is parallel with x-axis";
    if Degree(Ti,U) gt 1 then
      Append(~Tis,Ti);
    end if;
  end for;
  if #Tis gt 0 then
    return composite(Tis,U);            //SCHLECHT !! Besser KKT oder Kummer-T.
  else
    return U;
  end if;

end intrinsic;



//////////////////////////ordering field extension of local fields totally ramified over unramified/////////////////
intrinsic local_to_padic(L::FldPad,K::FldPad)->.{arrange the order of towers ob subfield so that totally ramified is over unramified}
           poly := DefiningPolynomial(L,K);
           loc := LocalField(K,poly);
           R,m := RamifiedRepresentation(loc);
           return R;
end intrinsic;



intrinsic arrange_field(L::FldPad,K::FldPad)->.{}
          if RamificationDegree(L,BaseField(L)) gt 1 and InertiaDegree(BaseField(L),K) eq Degree(BaseField(L),K) then
        return L;
//      end if;

        elif RamificationDegree(L,BaseField(L)) gt 1 and RamificationDegree(BaseField(BaseField(L)), K) gt 1  then
          return local_to_padic(L,K);
         //end if;  
        elif InertiaDegree(L,BaseField(L)) gt 1 and RamificationDegree(BaseField(L), K) gt 1 then
           return local_to_padic(L,K);
        else return L;
      end if;
end intrinsic;
intrinsic arrange_field(L::RngPad,K::RngPad)->.{}
          if RamificationDegree(L,BaseRing(L)) gt 1 and InertiaDegree(BaseRing(L),K) eq Degree(BaseRing(L),K) then
        return L;
//      end if;

        elif RamificationDegree(L,BaseRing(L)) gt 1 and RamificationDegree(BaseRing(BaseRing(L)), K) gt 1  then
          return local_to_padic(L,K);
         //end if;  
        elif InertiaDegree(L,BaseRing(L)) gt 1 and RamificationDegree(BaseRing(L), K) gt 1 then
           return local_to_padic(L,K);
        else return L;
      end if;
end intrinsic;




intrinsic factorise_easy(f,T,s)->.{easy way of factorising over small degree extension and then go for check up to extension over L1}

        K := CoefficientRing(f);
        assert #s eq 2;
        g := Expand(s[2]);
        _,_,E := Factorization(Polynomial(K,g):Extensions:=true);
        T1 := E[1]`Extension;
        _,_,E := Factorization(Polynomial(T,g):Extensions:=true);
        F := E[1]`Extension;
        f_d := InertiaDegree(F,PrimeField(F));
        fac := Factorisation( Polynomial(T1,f));
        A :=[];
        for x in fac do
            Append(~A,x[1]);
        end for;
        T := arrange_field(T,PrimeField(T));
        e_pol := DefiningPolynomial(T,BaseField(T));
        _,_,E := Factorization(Polynomial(T1,e_pol):Extensions:=true);
        T2 := E[1]`Extension;
        AA:=[];
        for x in A do
           fac_x := Factorisation(Polynomial(T2, x));
           for y in fac_x do
                Append(~AA,y[1]);
           end for;
        end for;
	F_equiv := ext<T2|f_d>;
assert AbsoluteDegree(F_equiv) eq AbsoluteDegree(F);

return [Polynomial(F_equiv, x) : x in AA], F_equiv;
end intrinsic;



intrinsic factorise_easy1(f,T,s)->.{easy way of factorising over small degree extension and then go for check up to extension over L1}
	
	T := arrange_field(T,PrimeField(T));
        K := CoefficientRing(f);
        assert #s eq 2;
        g := Expand(s[2]);
  //      _,_,E := Factorization(Polynomial(K,g):Extensions:=true);
//        T1 := E[1]`Extension;
        //_,_,E := Factorization(Polynomial(T,g):Extensions:=true);
        //F := E[1]`Extension;
        f_d := InertiaDegree(T,PrimeField(T));
 //       T := arrange_field(T,PrimeField(T));
        if AbsoluteDegree(T) eq InertiaDegree(T) then
           Unr := T;
	elif  AbsoluteDegree( BaseField(T)) eq  InertiaDegree(BaseField(T)) then
            Unr := BaseField(T);
         else Unr := ext<K|f_d>;// have to check here;
end if;
       // Unr := ext<K|f_d>;
        _,_,E := Factorization(Polynomial(Unr,g):Extensions:=true);
        T1 := E[1]`Extension;

       fac,_,fac_ext := Factorisation( Polynomial(T1,f): Extensions :=true);
       assert #fac ge 2;
       A :=[];
	B :=[];
        for i in [1..#fac] do
            Append(~A,fac[i,1]);
	    Append(~B,fac_ext[i]);
        end for;
        //T := arrange_field(T,PrimeField(T));
        e_pol := DefiningPolynomial(T,BaseField(T));
      assert RamificationDegree(T,PrimeField(T)) eq Degree(e_pol);
      _,_,E := Factorization(Polynomial(T1,e_pol):Extensions:=true);
	assert #E eq 1;//"irreducibility";
	T2 := E[1]`Extension;
       // A:=[Polynomial(T2,x[1]): x in A];
  E_pol,_,E := Factorisation(Polynomial(T2, A[1]): Extensions :=true);
  E_pol := [x[1]: x in E_pol];       
//        T2 := RamifiedRepresentation(LocalField(T1, Polynomial(T1,e_pol)));
     /*   AA:=[];
	BB := [];
        for x in A do
           fac_pol,_,fac_x := Factorisation(Polynomial(T2, x): Extensions :=true);
           for i in [1..#fac_x] do
                Append(~AA,fac_x[i]);
		Append(~BB, fac_pol[i]);
           end for;
        end for;*/
       // F_equiv := T2;  //ext<T2|f_d>;
//assert AbsoluteDegree(F_equiv) eq AbsoluteDegree(F);
return E_pol,T2, E;
//return [Polynomial(T2, x) : x in AA],T2;
end intrinsic;

/*//"the following does not work";
intrinsic gcd_facotrs(A::SeqEnum,f::RngUPolElt)->RngUPolElt
{ A consists of factors of f and it computes the GCD of (&*(A), f), alternatively  }
      for a in A do
	   f:=f div a;
       end for;
return f;
end intrinsic;
*/



intrinsic orbit_factor(psi::Map, f::RngUPolElt, g::RngUPolElt)-> SeqEnum,SeqEnum
{For g one factor of f, find all the orbit elements from subgroup of Domain(psi) so that psi(g)`s are different and hence we find all factors of (f), modified versio of ,,orbit_factor_G,,}
        G := Domain(psi);
        gen := [x : x in Generators(G)];
        Act := [];  //"create different orbit elements of G";
        for x in gen do
            for e in Orbit(G,x) do
                Append(~Act, e);
             end for;
         end for;
        F := CoefficientRing(g);
        n := Degree(f) div Degree(g);
        R<y> := PolynomialRing(F);
        h := hom<G->R| gamma :-> Polynomial([x@ psi(gamma) : x in Eltseq(g) ]) >;
        pols := [h(Id(G))];
        orbit := [Id(G)];
        for x in Act do
           delta := h(x);
           if delta notin pols then
              Append(~pols, delta);
              Append(~orbit, x);
           end if;
        if #pols ge n then
             gcd := GCD(&*(pols),R!f);
            if  Degree( gcd) eq Degree(f) then
                break x;
           // else return orbit_factor_G( psi, f, g);
            end if;
        end if;
        end for;
	set_minus := [x: x in G | not x in Act];
        if not  Degree(&*pols) eq Degree(f) then
	   for x in set_minus do
	       delta := h(x);
               if delta notin pols then
                  Append(~pols, delta);
                   Append(~orbit, x);
               end if;
	       if #pols ge n then
                  gcd := GCD(&*(pols),R!f);  //"This is expensive and it is just to claim that we have got all factors";
                  if  Degree( gcd) eq Degree(f) then
                      break x;
           // else return orbit_factor_G( psi, f, g);
                  end if;
              end if;
	  end for;
         // return orbit_factor_G(psi,f, g);
        end if;
return orbit, pols;
end intrinsic;


intrinsic orbit_factor_check(psi::Map, f::RngUPolElt, g::RngUPolElt)-> SeqEnum,SeqEnum
{For g one factor of f, find all the orbit elements from subgroup of Domain(psi) so that psi(g)`s are different and hence we find all factors of (f), without using GCD}
        G := Domain(psi);
        gen := [x : x in Generators(G)];
	Qp := CoefficientRing(f);
	n_0 := Precision(Qp);
	 F := CoefficientRing(g);
	FF:= ChangePrecision(F,n_0);
        Act := [];  //"create different orbit elements of G";
        for x in gen do
            for e in Orbit(G,x) do
                Append(~Act, e);
             end for;
         end for;
        //F := CoefficientRing(g);
        n := Degree(f) div Degree(g);
        R<y> := PolynomialRing(F);
        h := hom<G->R| gamma :-> Polynomial([x@ psi(gamma) : x in Eltseq(g) ]) >;
        pols := [h(Id(G))];
        pols_prec := [Polynomial(FF,pols[1])];
	orbit := [Id(G)];
        for x in Act do
           delta := h(x);
          // if delta notin pols then
           if Polynomial(FF,delta) notin pols_prec then   
	      Append(~pols, delta);
	      Append(~pols_prec, Polynomial(FF,delta));
              Append(~orbit, x);
           end if;
        if #pols ge n then
             gcd := GCD(&*(pols),R!f);
            if  Degree( gcd) eq Degree(f) then
                break x;
           // else return orbit_factor_G( psi, f, g);
            end if;
        end if;
        end for;
	
        set_minus := [x: x in G | not x in Act];
        if not  Degree(&*pols_prec) eq Degree(f) then
           for x in set_minus do
               delta := h(x);
               //if delta notin pols then
               if Polynomial(FF,delta) notin pols_prec then
		  Append(~pols, delta);
		  Append(~pols_prec, Polynomial(FF,delta));
                  Append(~orbit, x);
               end if;
               if #pols ge n then
                  gcd := GCD(&*(pols),R!f);  //"This is expensive and it is just to claim that we have got all factors";
                  if  Degree( gcd) eq Degree(f) then //"Alternatively,&and([Valuation(a) ge n_0: a in Coefficients(&*(pols)-R!f)]) eq true ";
                      break x;
           // else return orbit_factor_G( psi, f, g);
                  end if;
              end if;
          end for;
         // return orbit_factor_G(psi,f, g);
        end if;
return orbit, pols, pols_prec;
end intrinsic;



intrinsic orbit_factor_G(psi::Map, f::RngUPolElt, g::RngUPolElt)-> SeqEnum,SeqEnum
{For g one factor of f, find all the orbit elements from Domain(psi) so that psi(g)`s are different and hence we find all factors of (f)}
	G := Domain(psi);
	F := CoefficientRing(g);
	n := Degree(f) div Degree(g); 
	R<y> := PolynomialRing(F);
	h := hom<G->R| gamma :-> Polynomial([x@ psi(gamma) : x in Eltseq(g) ]) >;
	pols := [h(Id(G))];
	orbit := [Id(G)];
	for x in G do 
	   delta := h(x);
	   if delta notin pols then
	      Append(~pols, delta);
	      Append(~orbit, x);
	   end if;
	if #pols ge n and Degree(GCD(&*(pols), R!f)) eq Degree(f) then 
	  break x;
	end if;
	end for;
return orbit, pols;
end intrinsic;

intrinsic TwoSegmentsGaloisGroup_fast(f::RngUPolElt) -> GrpFP
  {Zu einem Eisensteinpolynom mit zwei Segmenten im VP wird die Galoisgruppe als endlich präsentierte Gruppe berechnet.}

 Qp:=CoefficientRing(f);
 p:=Prime(Qp);
  n:=Degree(f);
  e0:=n div p^Valuation(n,p);
  //1;
  T:=pReduction_e0(f);
  T := arrange_field(T,PrimeField(T));
  traegheit:=InertiaDegree(T,Qp);
  //1;
  if e0 eq 1 then
        s:=RamPolySubfields(f);
       assert Type(s) eq SeqEnum;//"this assertions fails if the precision of polynomial is not good enough";

       if Type(s) eq MonStgElt then 
          return "increase the precision of polynomial ring";
	end if;  
       g:=Expand(s[2]);
      _,_,E:=Factorization(Polynomial(T,g):Extensions:=true);
    F:=E[1]`Extension;
else          //gemischter Grad, TK zum 2. Segment ist in T enthalten!
F:=pReduction(f);
 end if;

        tt:=OmTree(Polynomial(Integers(F),f));//"prof. Pauli's code";
        vprint Segment, 1: "Computation time of factorisation of f";
       // vtime group, 1: E_pol:=[OmLift(t, Precision(Qp)): t in tt];//"only one factor is enought to do computation";
 	vtime group, 1: E_pol:=[OmLift(tt[1], Precision(Qp))];
	//  vtime group, 1: E_pol:=[OmLift(t, RamificationDegree(F,Qp)): t in tt];
 	E_pol:= [Polynomial(F,x): x in E_pol];
 prec:=UnitPrecision(Integers(F));
/*
 e_F := RamificationDegree(F,Qp);
 prec_F := e_F + prec+1 + 3+10;// "for extra if frovb equation loose the precision";
  
  */

//vprint Segment, 1: "Computation time of lfc";
//vtime group, 1: gamma,psi:= LFC_Using_Hilbert(F,Qp,prec+1);

//gamma,psi:= CLocalFundamentalClassSerre_check(F,Qp,prec+1); //es wird der Kozykel für Operation von links berechnet !!!!!!!!
  //3;
  //4;
 //4;
 /*      d1 := Degree(E_pol[1]);
        assert IsDivisibleBy(Degree(f), d1);
        act_d := Degree(f) div d1;*/
       // act_G := [x: x in G| Order(x) le act_d ]; // "this is not optimal set";
       //vtime group, 1: E_pol:= [x[1]: x in Factorisation(Polynomial(F,f))];
	
	 _,_,E:=Factorization(Polynomial(F,E_pol[1]):Extensions:=true);
         assert #E eq 1;
         Ext_1 := E[1]`Extension;
//Ext_1:= RamifiedRepresentation(LocalField(F,E_pol[1]));
       /* e_F := RamificationDegree(F,Qp);
        prec_F := e_F + prec+1 + 3+10;// "for extra if frovb equation loose the precision";

        if prec_F le Precision(F) then
           F := ChangePrecision(F,prec_F);
        end if;*/
/*
     disc_v := Valuation(Discriminant(F,K)) div InertiaDegree(F,K);//"this works!!!";
     prec_req := disc_v + 5 + prec+Degree(F);
     FF := ChangePrecision(F,prec_req);
     vprint Segment, 1: "Computation time of lfc";
     vtime group, 1: gamma,psi:= LFC_Using_Hilbert(F,PrimeField(F),prec+1);
     G:=Domain(gamma)[1];


*/
        vprint Segment, 1: "Computation time of lfc";
	vtime group, 1: gamma,psi:= LFC_Using_Hilbert(F,PrimeField(F),prec+1);
	G:=Domain(gamma)[1];

	 gamma:=map<Domain(gamma) -> Codomain(gamma) | x :-> gamma(x[2]^(-1),x[1]^(-1))>; //"Rechts-Kozykel"
	  //3;
         prec_norm := (prec+2)*Degree(Ext_1);
         vtime group, 1: orbit, pols := orbit_factor_check(psi, f, E_pol[1]);

         U,mU:=UnitGroup(F:Prec:=prec); //U;
	 m:=#Generators(U);

	vprint Segment, 1: "Computation time of norm group and its intersection of all NormGroup of factorisation's extension";
//vtime group, 1: norm_grp,map := NormGroup(Ext_1, mU);
        vtime group, 1: norm_grp,map := NormGroup_fast(Ext_1,mU, prec_norm);
        assert #quo<U|norm_grp> eq Degree(Ext_1,F); // " being abelian extensions";
	vtime group, 1: S:=&meet[ sub<U| [x@ map @mU@psi(g)@@mU: x in Generators(norm_grp)]>: g in orbit];


//vprint Segment, 1: "Computation time of intersection NormGroup of factorisation's extension";
//vtime group, 1: S:=&meet[ NormGroup(E[i]`Extension,mU) : i in [1..#E] ]; //S;
  US,mUS:=quo<U|S>;
  d:=#Generators(US);

  if S subset sub<U|[U.i : i in [1..m-1]], p*U.m> then  //ClassF zu S enthält die unverzw. Erw.
    traegheit:=traegheit*p;
  end if;

  //Mats:=[ Matrix([ Eltseq(mUS(Inverse(mU)(psi(G.j)(mU(Inverse(mUS)(US.i)))))) : i in [1..d] ])
  //       : j in [1..#Generators(G)+1] ];

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
  // Ext_F := Extension(CM,IdentifyTwoCocycle(CM,Gamma));
  // Ext_M := Group<S| S^I_d =1>;

 comp := rec<GaloisGrp | p:=p, gamma := gamma, E := E_pol, S:= S>;
 set := [* gamma, CM, mU, mUS, pols, F, psi, norm_grp, comp *];
 return //[ IdentifyGroup(GG) : GG in DistinctExtensions(CM) ];
         Extension(CM,IdentifyTwoCocycle(CM,Gamma)), traegheit, set, comp;

end intrinsic;
