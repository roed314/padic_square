//- MySplittingField: Zerfaellungskoerper ueber den Trivialansatz
//- DDGaloisGroup:    Galoisgruppe        -          "          - (automorphismGroup von Ruben)
//- pReduction: f in K[x] eisenstein mit p teilt Grad(f)-> Erweiterung T/K, so dass N/T p-Erweiterung ist
//- NormGroupSF: f(x) in K[x] eisenstein, Grad(f)=p^m -> Normgruppe des Zerfaellungskoerpers N, Abschaetzung fuer [N:K]


//////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 5.3 
intrinsic MySplittingField(f::RngUPolElt) -> FldPad
{Zerfaellungskörper von f (Trivialansatz)}
	
	while f ne One(PolynomialRing(CoefficientRing(f))) do
		
		fak,_,ext:=Factorization(f:Extensions:=true);
		1;
		lf:=[l : l in fak | Degree(l[1]) eq 1];
		for i in [1..#lf] do
			f:=f div lf[i][1];
		end for;
		_,max:=Maximum([Degree(a[1]) : a in fak]);
		L:=ext[max]`Extension;
		
		//if exists(n){n : n in [1..#fak] | Degree(fak[n][1]) gt 1} then
		//	L:=ext[n]`Extension;
		//end if;
		
		f:=Polynomial(L,f);
		//pr:=SuggestedPrecision(f);pr;
		//ChangePrecision(~L,30);
		
	end while;
	
	return L;

end intrinsic;


intrinsic DDGaloisGroup(f::RngUPolElt) -> Grp, Map, PowMapAut
{Galoisgruppe von f als Permutationsgruppe (Trivialansatz)}
	//return AutomorphismGroup(MySplittingField(f),CoefficientRing(f));
	return AutomorphismGroup( Integers(MySplittingField(f)), Integers(CoefficientRing(f)) ); //siehe lfc.m
end intrinsic;



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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 5.1
intrinsic pReduction(f::RngUPolElt) -> FldPad
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
  vpInfos:=XFlatRamificationPolygon(f);	// [ ..., < h_i, e_i, p^{s_{i-1}}, A_i(y)>, ... ]
  if flaches_Segment then
    ell:=#vpInfos-1;	//letztes Segment separat behandeln
  else
    ell:=#vpInfos;
  end if;  
  
  b:=[];
  v:=[];
  eff:=[];
  for i in [1..ell] do		//letztes Segment separat behandeln!
    _,bi,_:=XGCD(vpInfos[i][1],vpInfos[i][2]);
    Append(~b,bi);
    vi:=bi*(n div vpInfos[i][3])+n+1;
    Append(~v,vi);
    f1:=LCM([Degree(s[1]) : s in Factorization(vpInfos[i][4])]); //ass. Trägheit
    f0:=degree_zeta_e(e0*vpInfos[i][2],q);
    Append(~eff,LCM(f0,f1));
  end for;
  
  //if flaches_Segment then
  //  Append(~eff,degree_zeta_e(e0,q));
  //end if;  
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
    Ti:=ext<U|x^(e0*vpInfos[i][2])-((-1)^v[i]*eps^(b[i]*n)*Coefficient(f,0))>;
    if Degree(Ti,U) gt 1 then
      Append(~Tis,Ti);
    end if;
  end for;
  
  if #Tis gt 0 then
    return composite(Tis,U);		//SCHLECHT !! Besser KKT oder Kummer-T.
  else 
    return U;
  end if;
  
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

intrinsic CopyNormGroup(U::GrpAb, iso::Map, N::GrpAb, aut::Map) -> GrpAb
	{}
	gens:=[];
	for n in Generators(N) do
		m:=Inverse(iso)(aut(iso(n)));
		Append(~gens,m);
	end for;

	return sub<U|gens>;
end intrinsic;


intrinsic Intersec(list::SeqEnum) -> GrpAb
	{}
	S:=list[1] meet list[2];
	for i in [3..#list] do
		S:=S meet list[i];
	end for;
	
	return S;
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ACHTUNG: - Berücksichtigt nur abelsche Schnitte der Teilbäume! -> Nur obere Schranke für Grad des Zerf`kps
//Rekursive Funktion fuer NormGroupSF()
intrinsic RecSF(f::RngUPolElt) -> RngIntElt, GrpAb, GrpAb, Map
  {}
  K:=CoefficientRing(f);
  V,mV:=UnitGroup(K);	//Präzision angeben
  d:=Degree(f);
  s:=RamPolySubfields(f);
  
  if #s eq 1 then	//1 segment
    L:=ext<K|f>;
    return Degree(L,K),NormGroup(L,mV),V,mV; 
  end if;
  
  for i in [2..#s] do
    L:=ext<K|Expand(s[i])>;
    fac,_,ext:=Factorization(Polynomial(L,f):Extensions:=true);
    deg:=[ Degree(fac[j][1]) : j in [1..#fac] ];deg;	//assertion ??
    
    if SequenceToSet(deg) eq { d div Degree(L,K) } then //L/K galoissch (!!!!!?)
      //U,mU:=UnitGroup(L);
      if i eq 2 then		//ganz oben im baum
	1;
	U,mU:=UnitGroup(L);	//Präzision ?
	1;
	N:=Intersec([NormGroup(e`Extension,mU) : e in ext]);// U/N entspr. ab.Q.(= komplette gruppe) über L
	//U/N;
	n:=Degree(L,K)*#(U/N);
	N:=sub<V| [Inverse(mV)(Norm(mU(u))):u in Generators(N)]>; // V/N entspr. ab.Q. über K
	//V/N;
	return n,N,V,mV;
      end if;
      
      n,N,U,mU:=RecSF(fac[1][1]);	//fac[1][1] wirklich ausreichend? annahme: die galoisgruppen von ext[i]/K sind isomorph !!
      a:=#(U/N); //a;//grad des abQ für einen ast
      Aut:=Automorphisms(L,K); 			//assertion: #Aut eq Degree(L,K)
      N:=Intersec([ CopyNormGroup(U,mU,N,Aut[j]) : j in [1..#Aut] ]); //(U/N);
      as:=(a^Degree(L,K)) div #(U/N);//faktor, der durch den schnitt der abQ verloren geht
      N:=sub<V| [Inverse(mV)(Norm(mU(u))):u in Generators(N)]>; //s.o, V/N entspr. ab.Q. der normalen hülle von L/K
      n:=(Degree(L,K)*n^Degree(L,K)) div as; 
      
      return n,N,V,mV;
    end if;  
  end for;  
      
end intrinsic;  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Nur für Grad-p^m-Eisensteinpolynome!!
//vgl. Algorithmus 5.5
intrinsic NormGroupSF(f::RngUPolElt) -> GrpAb, RngIntElt
  {f(x)\in K[x] eisenstein, Grad(f)=p^m -> Normgruppe des Zerfaellungskoerpers N, Abschaetzung fuer [N:K]}
  
  K:=CoefficientRing(f);
  V,mV:=UnitGroup(K); //Präzision!
  T:=pReduction(f);
  if RamificationIndex(T,K) gt 1 then
    f:=Eisenstein(Polynomial(T,f));
  else
    f:=Polynomial(T,f);
  end if;
  
  //Vorsicht: Bei Grad e_0*p^m kann T kopieren!!
  n,N,U,mU:=RecSF(f);	//über T ...
  n:=n*Degree(T,K);
  N:=sub<V| [Inverse(mV)(Norm(mU(u),K)):u in Generators(N)]>; //über K ...
  
  return V/N, n;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
intrinsic NewCFSplittingField(f::RngUPolElt) ->  FldPad
	{ACHTUNG: Noch nicht fertig!!}
	
	Kf:=CoefficientRing(f);
	n:=Degree(f);
	p:=Prime(Kf);
	e0:=n div p^Valuation(n,p);
	s:=Reverse(RamPolySubfields(f));
	K:=pReduction(f);
	
	if e0 eq 1 then
	  K:=ext<K|Eisenstein(Polynomial(K,Expand(s[1])))>;	//geschickter lösen !!
	end if;  

	for i in [2..#s] do
		//KX<x>:=PolynomialRing(K);
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);//E;
		//cfs:=[E[j]`Extension : j in [1..#E]];
		pr:=UnitPrecision(Integers(K));
		1;
		U,mU:=UnitGroup(K :Prec:=pr);		//pSelmerGroup ??
		2;
		//normgrps:=[];
		erz:=[];
		Append(~erz, E[1]`Extension);
		N:=NormGroup(E[1]`Extension,mU);
		for j in [2..#E] do			//ACHTUNG! Annahme: Zwei versch. "Kopien" sind entweder gleich oder disjunkt.
		  Nj:=NormGroup(E[j]`Extension,mU);	//F A L S C H !!!!!!!!!!!!!
		  if not (N subset Nj) then
		    Append(~erz, E[j]`Extension);
		    N:=N meet Nj;
		  end if;  
		end for;  
		
		//gi:=DefiningPolynomial(E[1]`Extension);       
		//pr:=UnitPrecision(Integers(K));
		//A:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		//if #A ne Degree(K,Kf) then
		//	error "Fehler in der Automorphismenberechnung!";
		//end if;	
		//if i ne 2 then ExtTest(U,iso,A,p); end if;
		//N:=NormGroup(Ki,iso);
		//normgrps:=[];
		//for aut in A do 				//ein EZS würde auch reichen, viele aut operieren trivial!
		//	M:=CopyNormGroup(U,iso,N,aut);
		//	Append(~normgrps,M);
		//end for;
		
		
		//NN:=Intersec(normgrps);
		//cfs:=ClassField(mU,NN);
		//3;
		//K:=CFComposite(cfs);
		//4;
	end for;
    
	return erz; //K;

end intrinsic;*/
