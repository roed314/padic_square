//RamPolyFactors: f(x) eisenstein -> Faktorisierung von f(x) zum Verzweigungspolygon
//RamPolySubfields: f(x) eisenstein -> Minimalpolynome der kanonischen Teilkörper zum Verzweigungspolygon
//CFSplittingField: Zu einem Eisensteinpolynom wird der Zerfaellungskoerper anhand des Teilkoerperturms zum VP berechnet.
//                  Es wird Klassenkörpertheorie zum "Kopieren" benutzt 
//CSplittingField: Wie CFSplittingField, allerdings ohne Klassenkörpertheorie.


intrinsic Expand(f::RngUPolElt) -> RngUPolElt
	{}
	L:=[ Expand(a) : a in Coefficients(f) ];
	return Polynomial(CoefficientRing(f),L);
end intrinsic;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 4.3
intrinsic RamPolyFactors(f::RngUPolElt) -> SeqEnum //, Map
	{f(x) eisenstein -> Faktorisierung von f(x) "zum Verzweigungspolygon"}

	//Lnew<a>:=LocalField(CoefficientRing(f),f);      //für RamPolySubfields(), Berechnung der Teilkörper in RngLocA leichter,
	//L<alpha>,iso:=RamifiedRepresentation(Lnew);     //aber sub<> funktioniert nicht immer...UND ALPHA IST NICHT NS VON f!!

	L<alpha>:=ext<CoefficientRing(f)|f>;

	P<x>:=PolynomialRing(L);
	rho:=Evaluate(f,alpha*x+alpha) div (x*alpha^Degree(f));	
	v:=Vertices(NewtonPolygon(rho));
	factors:=[];

	for i in Reverse([2..#v]) do
		N:=Degree(rho);
		m:=Integers() ! v[i-1][1];
		E:=Integers() ! (v[i][1]-v[i-1][1]);
		H:=Integers() ! (v[i-1][2]-v[i][2]);
		e:=E div GCD(E,H);
		h:=H div GCD(E,H);
		
		LL<beta>:=ext<L|x^e-alpha>;		
		R<y>:=PolynomialRing(LL);
		srho:=Evaluate(rho,y*beta^h) / beta^(N*h);//srho;                       //"flatten" the last segment
		c:=[ Coefficient(rho,j)*beta^(h*(j-N)) : j in [m..N] ];
		srho2:=Polynomial(LL,c);
		srho1:=y^m;
		//11;
		hl:=HenselLift(Polynomial(Integers(LL),srho),[srho1,srho2]);//hl;
		//11;
		rho1:=Polynomial(L,Evaluate(hl[1],y/beta^h)*(beta^h)^m);
		rho2:=Polynomial(L,Evaluate(hl[2],y/beta^h)*(beta^h)^(N-m));    //should be over L (unique factor cor. to last segment)
		Append(~factors,rho2);
		
		rho:=rho1;
	end for;

	factors:=Reverse(factors);//factors;
	factors:=[ Evaluate(rhoi,(x-alpha)/alpha)*alpha^Degree(rhoi) : rhoi in factors ];    //back to f...
	factors[1]:=factors[1]*(x-alpha);

	return factors,L; //, iso;

end intrinsic;


////////////////////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 4.4
intrinsic RamPolySubfields(f::RngUPolElt) -> SeqEnum  	
	{f(x) eisenstein -> Minimalpolynome der kanonischen Teilkörper zum Verzweigungspolygon}
	
	K:=CoefficientRing(f);
	fak:=RamPolyFactors(f);
	list:=[]; Append(~list,f);
	
	/*c:=Coefficient(fak[1],0);
	m:=MinimalPolynomial(c,K);
	if Degree(m) ne 
	Append(~list, m);*/

	c:=1;
	for i in [1..#fak-1] do 
		c:=c*Coefficient(fak[i],0);
		minpo:=MinimalPolynomial(c,K);
		if Degree(minpo) eq Degree(f) then		//MinimalPolynomial() berechnet manchmal nur char. Pol.!
			print "RamPolySubfields: Minimalpolynom musste korrigiert werden!";
			minpo:=Factorization(minpo)[1][1];      //häufig Präzisionsprobleme da nicht quadratfrei...
		end if;
		Append(~list, minpo);  
	end for;

	return list;
	
end intrinsic;

////////////////////////////////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 4.5
//FEHLERHAFT: Die lineare Algebra über K funktioniert nur selten ...
intrinsic RamPolyTower(f::RngUPolElt) -> FldPad
  {}
  K:=CoefficientRing(f);
  fac,L:=RamPolyFactors(Expand(f));
  if #fac eq 1 then
    return L;
  end if;
  h:=&*[fac[i] : i in [1..#fac-1]]; //h in L[x] ("Relativpolynom")
  Li:=ext<K|Expand(MinimalPolynomial(Coefficient(h,0),K))>; 
  //h in Li[x] interpretieren:
  M:=Matrix(K,[ Eltseq(Coefficient(h,0)^i) : i in [0..Degree(Li,K)-1] ]);
  v:=[ Vector(K,Eltseq(Coefficient(h,i))) : i in [0..Degree(h)] ];
  s:=[ Solution(M,v[i]) : i in [1..#v] ];			//ACHTUNG: Präzisionsprobleme!!
  P<x>:=PolynomialRing(Li);
  h:=P![Li!Eltseq(s[i]) : i in [1..#s]]; //h in Li[x]
  return RamPolyTower(h);
end intrinsic;




///////////aus class_field.m////////////////////////////////////////////////////
function eval_automorphism(x, gens, pos)
    if pos eq #gens then
        return x;
    end if;

    R := Parent(x);
    P := PolynomialRing(R);
    return Evaluate(
        P ! [eval_automorphism(i, gens, pos + 1) : i in Eltseq(x)], gens[pos]
    );
end function;


function eval_field_automorphism(x, m)
    F := Parent(x);
    R := RingOfIntegers(F);
    pi := UniformizingElement(F);

    xv := Valuation(x);
    xu := ShiftValuation(x, -xv);
    return (F ! m(UniformizingElement(R)))^xv * F ! m(R ! xu);
end function;

function construct_field_map(L, m)
    res := map<L -> L | x :-> eval_field_automorphism(x, m)>;
    res`LocalRingMap := m;
    return res;
end function;

function auto_eq(m1, m2)
    return m1`LocalAutoGenerators eq m2`LocalAutoGenerators;
end function;

function auto_mult(m1, m2)
    gens1 := m1`LocalAutoGenerators;
    gens2 := m2`LocalAutoGenerators;
    gens := [];
    for i in [1..#gens1] do
        Append(~gens, eval_automorphism(gens1[i], gens2, 1));
    end for;

    L := Universe(gens1);
    m := map<L -> L | x :-> eval_automorphism(x, gens, 1)>;
    m`LocalAutoGenerators := gens;
    return m;
end function;
/////////////////////////////////////////////////////////////////////////////////


intrinsic ConstructFieldMap(L::FldPad, m::Map) -> Map
	{}
	return construct_field_map(L,m);
end intrinsic;


//intrinsic Continuations(m::Map, L::FldPad) -> Map
//{Let m: R -> R be an automorphism of a p-adic field.
//Let L be an extension of R.  Return the continuations of m to L.}
//    maps := Continuations(m,RingOfIntegers(L));
//    return [construct_field_map(L, m) : m in maps];
//end intrinsic;


function auto_eq_field(m1, m2)
    return m1`LocalRingMap`LocalAutoGenerators eq m2`LocalRingMap`LocalAutoGenerators;
end function;

function auto_mult_field(m1, m2)
    m:=auto_mult(m1`LocalRingMap,m2`LocalRingMap);
    return construct_field_map(Domain(m1),m);
end function;



intrinsic AutToGroup(A::SeqEnum) -> GrpFP,Map
	{Fehlerhaft bei größeren Gruppen. GenericGroup?}
	return GenericGroup(A : Mult := auto_mult_field, Eq := auto_eq_field, Id := A[1]);
end intrinsic;
    


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
	


intrinsic CFComposite(list::SeqEnum) -> FldPad					//einfacher in RngLocA!!
	{Kompositum der von ClassField() ausgegebenen zyklischen Erweiterungen.}
	N:=list[1];
	for i in [2..#list] do
		_,_,E:=Factorization(Polynomial(N,DefiningPolynomial(list[i])):Extensions:=true);
		N:=E[1]`Extension;
	end for;
	return N;
end intrinsic;



intrinsic CFComposite(L::FldPad) -> FldPad
	{ClassField() gibt nicht immer eine Liste aus!}
	return L;
end intrinsic;


function Fe(R)
  p := Prime(R);
  e := RamificationIndex(R,PrimeRing(R));
  F := [];
  for i := 1 to Floor(p*e/(p-1)) do 
    if i mod p ne 0 then 
      Append(~F,i);
    end if;
  end for;
  return F;
end function;


intrinsic FLevels(R::RngPad) -> SeqEnum
	{}
	return Fe(R);
end intrinsic;
	


intrinsic UnitPrecision(R::RngPad) -> RngIntElt
	{}
	p:=Prime(R);
	e:=RamificationIndex(R,PrimeRing(R));
	F:=Fe(R);

	if not HasPRoot(R) then 	
		return Maximum(F)+1;
	else
		mu0:=Valuation(e,p)+1;
		e0:= e div (p^(mu0-1)*(p-1));            // es gilt (p-1)|e!!
		return Maximum(F cat [p^mu0*e0])+1;
	end if;
end intrinsic;



intrinsic UnitTest(K::FldPad, L::FldPad) -> GrpAb, GrpAb, Map, GrpAb, Map
	{}
	p:=Prime(K);
	U1,m1:=PrincipalUnitGroup(Integers(K): Prec:=UnitPrecision(Integers(K))); 
	U2,m2:=PrincipalUnitGroup(Integers(L): Prec:=UnitPrecision(Integers(L)));
	F,h:=quo<U2|p*U2>;
	gens:=[ h(Inverse(m2)(m1(U1.i))) : i in [1..#Generators(U1)] ];
	S:=sub<F|gens>;
	return S,F,h,U2,m2;
end intrinsic;




intrinsic UnitTest2(K::FldPad, L::FldPad) -> GrpAb, GrpAb, Map, GrpAb, Map
	{}
	p:=Prime(K);
	U1,m1:=UnitGroup(K: Prec:=UnitPrecision(Integers(K))); 
	U2,m2:=UnitGroup(L: Prec:=UnitPrecision(Integers(L)));
	F,h:=quo<U2|p*U2>;
	gens:=[ h(Inverse(m2)(m1(U1.i))) : i in [1..#Generators(U1)] ];
	S:=sub<F|gens>;
	return S,F,h,U2,m2;
end intrinsic;




intrinsic ExtTest(U::GrpAb,mU::Map,A::SeqEnum,p::RngIntElt) -> RngIntElt
	{}
	QU,mQU:=quo<U|p*U>;
	n:=#Generators(QU);
	phi:=Inverse(mU)*mQU;
	mats:=[];
	for i in [2..#A] do
		op:=[ Eltseq( phi(A[i](Inverse(phi)(QU.j))) ) : j in [1..n] ];
		Append(~mats, GL(n,p)!Matrix(op));
	end for;
	Op:=MatrixGroup<n,GF(p)|mats>;
	M:=GModule(Op);
	CM:=CohomologyModule(Op,M);

	//return Op;
	return #DistinctExtensions(CM);                  
	//return CohomologyGroup(CM,2);
end intrinsic;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////77
//vgl. Algorithmus 5.4
intrinsic CFSplittingField(f::RngUPolElt) ->  FldPad
	{Zu einem Eisensteinpolynom wir der Zerfaellungskoerper anhand des Teilkoerperturms zum VP bestimmt. 
	 Dabei wird die Normgruppe kopiert und der Zerfällungskörper mit ClassFields() berechnet.}
	
	Kf:=CoefficientRing(f);
	p:=Prime(Kf);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;
	L:=K;

	for i in [2..#s] do
		//KX<x>:=PolynomialRing(K);
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        
		pr:=UnitPrecision(Integers(K));
		U,iso:=UnitGroup(K : Prec:=pr);
		A:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	
		//if i ne 2 then ExtTest(U,iso,A,p); end if;
		N:=NormGroup(Ki,iso);
		normgrps:=[];
		for aut in A do 				
			M:=CopyNormGroup(U,iso,N,aut);
			Append(~normgrps,M);
		end for;
		NN:=Intersec(normgrps);
		cfs:=ClassField(iso,NN);
		K:=CFComposite(cfs);
	end for;
    
	return K; //ExtTest(U,iso,A,p);

end intrinsic;


//////////////////////////////////////////////////////////////////////////////////////////////////
//vgl. Algorithmus 5.4
//Nur einfaches Kopieren der Polynome ohne Klassenkörpertheorie!
intrinsic CSplittingField(f::RngUPolElt) ->  FldPad
	{Zu einem Eisensteinpolynom wir der Zerfaellungskoerper anhand des Teilkoerperturms zum VP bestimmt.}
	
	Kf:=CoefficientRing(f);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;

	for i in [2..#s] do
		KX<x>:=PolynomialRing(K);
		1;
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		1;
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        //immer richtiger Faktor??
		2;
		A:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	
		2;
		3;
		L:=Ki;
		for aut in Remove(A,1) do 				
			hi:=CopyPolynomial(gi,aut);
			hi:=Polynomial(L,hi);
			if not HasRoot(hi) then                               //Bedingung korrekt??
				4;
			        _,_,EE:=Factorization(hi:Extensions:=true);
		                L:=EE[1]`Extension;
				4;
			end if;
		end for;
		3;
	        K:=L;	
	end for;
    
	return K;

end intrinsic;



intrinsic CFGaloisGroup(f::RngUPolElt) ->  FldPad, SeqEnum
	{}
	
	Kf:=CoefficientRing(f);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;

	A:=Automorphisms(K,Kf);                       
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	

	for i in [2..#s] do
		//KX<x>:=PolynomialRing(K);
		1;
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);  //bei 2 Segmenten ohne Fact?
		1;
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        //immer richtiger Faktor??
		2;
		U,iso:=UnitGroup(K: Prec:=UnitPrecision(Integers(K)));
		2;
		4;
		N:=NormGroup(Ki,iso);
		4;
		normgrps:=[];
		for aut in A do 				//ein EZS würde auch reichen!!
			M:=CopyNormGroup(U,iso,N,aut);
			Append(~normgrps,M);
		end for;
		NN:=Intersec(normgrps);
		5;
		cfs:=ClassField(iso,NN);
		5;
		6;	
		L:=CFComposite(cfs);
		6;
		7;
		AA:=[map<Integers(K)->Integers(L)|x:->aut(x)> : aut in A];
		cont:=[];
		for i in [1..#AA] do
			AA[i]`LocalAutoGenerators:=ChangeUniverse(A[i]`LocalRingMap`LocalAutoGenerators, Integers(L));
			cont cat:=[ConstructFieldMap(L,m) : m in Continuations(AA[i],Integers(L))];
		end for;
		7;
		A:=cont;
		K:=L;
	end for;
    
	return K,A;

end intrinsic;



intrinsic Q2SplittingField(f::RngUPolElt) ->  FldPad, SeqEnum//Grp, Map, PowMapAut, FldPad
	{Testversion fuer Erweiterungen von Q2, bei denen nur der kleinste Teilkörper nicht Grad 2 hat
	 - nur korrekt, wenn der Zerfällungskörper nur durch "Kopieren" bestimmt werden kann.}
	
	Kf:=CoefficientRing(f);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;

	for i in [2..#s] do
		KX<x>:=PolynomialRing(K);
		1;
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		1;
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        //immer richtiger Faktor??
		gi:=Evaluate(gi,x-(Coefficient(gi,1)/2));      //i.A. nicht eisenstein!
		
		pr:=UnitPrecision(Integers(K));
		pr;
		U,iso:=UnitGroup(K: Prec:=pr);
		QU,nhom:=quo<U|[2*u:u in Generators(U)]>;
		phi:=Inverse(iso)*nhom;                        //phi nicht injektiv!?
		2;
		A:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	
		2;
		L:=Ki;
		test:=[];
		Append(~test, phi(Coefficient(gi,0)));

		for aut in A do 
			a:=phi(aut(Coefficient(gi,0)));
			if a notin test then
				3;				                                              
				_,_,EE:=Factorization(Polynomial(L,x^2+Inverse(phi)(a)):Extensions:=true);
				3;
				//L:=ext<L|Expand(Polynomial(L,x^2-Inverse(phi)(a)))>; 
				L:=EE[1]`Extension;
				Append(~test, a);
			end if;
		end for;
		//L;
		
		if IsTotallyRamified(L) then
			K:=ext<Kf|DefiningPolynomial(L,Kf)>;         //besser für Aut-Berechnung??
		else
			K:=L;					     //falls unverzweigtes hinzu gekommen ist!!!
		end if;
	end for;
    
	return K; //Automorphisms(K,Kf); 

end intrinsic;




intrinsic GStrich(G::GrpFP, p::RngIntElt) -> GrpFP, Map
	{}
	pot:=[g^p : g in Generators(G)];
	com:=[];
	for i in [1..#Generators(G)-1] do
		for j in [i+1..#Generators(G)] do
			Append(~com, (G.i,G.j));
		end for;
	end for;	
	return sub<G|pot cat com>;     //oder ncl<>??
end intrinsic;



intrinsic Gab(G::GrpFP) -> GrpFP, Map
	{}
	com:=[];
	for i in [1..#Generators(G)-1] do
		for j in [i+1..#Generators(G)] do
			Append(~com, (G.i,G.j));
		end for;
	end for;	
	return ncl<G|com>;     //ncl<> klappt häufig nicht (no "closed coset table")
end intrinsic;
	


///////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------AB HIER ALTES UND TEILWEISE FEHLERHAFTES ZEUG-----------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////


//Bei SafarevicGG, FPGG und TSGG (?) werden zur Identifikation von K*/K*^2 mit R/R'R^2 keine Modulisos benutzt. Es wird kein Kozykel
//und auch kein Lenstra-Diagramm berücksichtigt!!



intrinsic SafarevicGaloisGroup(f::RngUPolElt) ->  GrpFP
	{Berechnet die Galoisgruppe als FP-Gruppe nach dem Ansatz von Safarevic. Der Zerfkp wird nicht vollst. konstruiert.
         Klappt z.Z. nur bei kleinstem Teilkörper vom Grad p und \zeta_p \notin Grundkörper und wenn die Galoisgruppe eine p-Gruppe ist.
  	 ???Wirklich korrekt bei mehr als 2 Segmenten???}
	
	Kf:=CoefficientRing(f);
	p:=Prime(Kf);
	n0:=Degree(Kf);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;
	if p eq 2 then
		F:=FreeGroup(3);
		Falt:=sub<F|F.1,F.2,F.3^2,(F.1,F.2),(F.1,F.3),(F.2,F.3)>;
	else
		F:=FreeGroup(n0+1);
		com:=[];
		for i in [1..n0] do
			for j in [i+1..n0+1] do
				Append(~com, (F.i,F.j));
			end for;
		end for;	
		Falt:=sub<F|[F.1^p] cat [F.i : i in [2..n0+1]] cat com>;                    //-bei größeren Segmenten C_p x...x C_p
	end if;

	for i in [2..#s] do		
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);
		pr:=UnitPrecision(Integers(K));
		//pr;
		U,iso:=UnitGroup(K : Prec:=pr);
		UF,mapUF:=quo<U|p*U>;
		n:=#Generators(UF);
		n;

		Aut:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		if #Aut ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	
		
		N:=NormGroup(Ki,iso);
	
		normgrps:=[];
		for aut in Aut do 				//ein EZS würde auch reichen, viele aut operieren trivial!
			M:=CopyNormGroup(U,iso,N,aut);
			Append(~normgrps,M);
		end for;
		NN:=Intersec(normgrps);

		phi:=Inverse(iso)*mapUF;                      
		OpM2:=[];
		for i in [2..#Aut] do
			op:=[ Eltseq( phi(Aut[i](Inverse(phi)(UF.j))) ) : j in [1..n] ];
			op;
			Append(~OpM2, GL(n,p)!Matrix(op));
		end for;
		OpM2:=MatrixGroup<n,GF(p)|OpM2>;                   //Operation auf der max el-ab Erw. aus der KKT (Automorphismen)
		OpM2;	

		Fneu:=GStrich(Falt,p);
		if p eq 2 then
			G:=quo<F|Fneu,F.1^2*F.2^4*(F.2,F.3)>;
		else
			G:=quo<F|Fneu>;
		end if;
		Order(G);
		A,mapA:=ncl<G|Generators(Falt)>;              //NT entspricht max el-ab Erweiterung
		Order(A);
		M1,mapM1:=GModule(G,A,p);

	/*	if p eq 2 then                                                   //////läuft nicht korrekt///////////
			Sub:=[ S : S in Submodules(M1) | Dimension(S) eq n ];
			Act:=[Image(GModuleAction(S)) : S in Sub ];
			bool:=false;
			k:=1;
			while not bool do 
				bool, B:=IsConjugate(GL(n,p),Act[k],OpM2);
				S1:=Sub[k]; k;
				k:=k+1;
			end while;
			//_,B:=IsConjugate(GL(n,p),Act[1],OpM2);
			//S1:=Sub[1];
			
			relsM2:=[ mapUF(U!NN.i) : i in [1..#Generators(NN)] ];
			relsM2:=[ Eltseq(a) : a in relsM2];
			//relsM1:=[ Vector(GF(p),a)*B : a in relsM2 ];    //Umrechnen der Relationen aus der Normgruppe
			relsS1:=[ Vector(GF(p),a)*B^-1 : a in relsM2 ];   //s.o.
		
			rels:=[];
			for v in relsS1 do
				r:=Id(G);
				for i in [1..Dimension(S1)] do 
					r:=r*mapA(Inverse(mapM1)(S1.i))^(Integers()!v[i]);   //mapA(Inverse(mapM1)(M1.1)) : M1 -> G
				end for;
				Append(~rels,r);
			end for;
			Q:=quo<M1|S1>;
			zrels:=[ mapA(Inverse(mapM1)(Q.i)) : i in [1..#Generators(Q)] ];
		
			G:=quo<G|rels cat zrels>;


		else */
			OpM1:=Image(GModuleAction(M1));		      //Operation auf dem NT aus der FP-Gruppe	
			OpM1;

			//_,B:=IsConjugate(GL(n,p),OpM2,OpM1);                   //M2^B=M1 Basiswechsel     //versch. Dimensionen bei p=2!!!!
			_,B:=IsConjugate(GL(n,p),OpM1,OpM2);	                 //in einigen Bsps schneller mit dieser Reihenfolge: M1^B=M2 -> M2^B^-1=M1

			relsM2:=[ mapUF(U!NN.i) : i in [1..#Generators(NN)] ];
			relsM2:=[ Eltseq(a) : a in relsM2];
			//relsM1:=[ Vector(GF(p),a)*B : a in relsM2 ];    //Umrechnen der Relationen aus der Normgruppe
			relsM1:=[ Vector(GF(p),a)*B^-1 : a in relsM2 ];   //s.o.
		
			rels:=[];
			for v in relsM1 do
				r:=Id(G);
				for i in [1..Dimension(M1)] do 
					r:=r*mapA(Inverse(mapM1)(M1.i))^(Integers()!v[i]);   //mapA(Inverse(mapM1)(M1.1)) : M1 -> G
				end for;
				Append(~rels,r);
			end for;
			
			
			G:=quo<G|rels>;
	      //end if;

		Falt:=Fneu;                                          //////RICHTIG????????//////////////////////

		if i ne #s then	
			cfs:=ClassField(iso,NN);           //im letzten Schritt keine explizite Konstr. des Körpers mehr nötig!
			K:=CFComposite(cfs);
		end if;
		
		
	end for;
    
	return G;

end intrinsic;





intrinsic TestPolynomial(K::FldPad, p::RngIntElt, n::RngIntElt) -> RngUPolElt
	{}
	P<x>:=PolynomialRing(K);
	OK:=Integers(K);
        //L:=[ [3,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,1] :
	//   a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26 in [0,3] ];
	//pols:=[];
	while true do
		u:=Random(OK);
		if Valuation(u) eq 0 then
			c:=[p*u];
			for i in [1..n-1] do
				//Append(~c,Random([0,3,9,27,6,12,15,3^4,3^5,3^6,3^7,3^8,3^9,3^10,3^11,2*3^4,2*3^5,2*3^6,2*3^7,2*3^8,2*3^9,2*3^10,2*3^11]));
				Append(~c, p*Random(OK));
			end for;
			Append(~c,1);
			f:=P ! c;
			f;
			r:=FlatRamificationPolygon(f);
			r;
			if #r eq 2 then							//wie allgemein für log_p(n)???
				//r;
				if r[1][1] eq 1 and r[2][1] eq 1 then                 //and r[3][1] eq 1 then
					return f;
				end if;
			end if;
		end if;
	end while;
end intrinsic;




intrinsic FPGaloisGroup(f::RngUPolElt, mat::AlgMatElt) -> GrpFP
	{}
	K:=CoefficientRing(f);
	p:=Prime(K);
	n0:=Degree(K);
	s:=Reverse(RamPolySubfields(f));
	L:=ext<K|Expand(s[1])>;

	if p eq 2 then
		FF:=FreeGroup(n0+2);
		F:=quo<FF|FF.1^2*FF.2^4*(FF.2,FF.3)>;
		//F:=quo<FF|FF.1^2*FF.3^4*(FF.3,FF.2)>;	
	else
		F:=FreeGroup(n0+1);
	end if;
	Falt:=F;

	for i in [1..#s] do
		//U,mU:=UnitGroup(K);               			//besser mit geringerer Präzision!!
		pr:=UnitPrecision(Integers(K));
		U,mU:=UnitGroup(K : Prec:=pr);
		UF,mUF:=quo<U|p*U>;
		n:=#Generators(UF);
		N,m:=NormGroup(L,mU);
		Fneu:=GStrich(Falt,p);
		/*if p eq 2 then
			Fneu:=ncl<F|Fneu,F.1^2*F.2^4*(F.2,F.3)>;	//zus. Rel. korrekt an dieser Stelle??
		end if;
		*/
		if i eq 1 then  
			relsN:=[ mUF(U!N.i) : i in [1..#Generators(N)] ];
			//relsN:=[ U!N.i : i in [1..#Generators(N)] ];	
			relsN:=[ Eltseq(a) : a in relsN];
			//relsN:=[[ y mod 2 : y in relsN[k] ] : k in [1..#relsN] ];
			//relsN:=[ Vector(GF(p),b)*Matrix(GF(p),[[0,0,1],[1,0,0],[0,1,0]]) : b in relsN ];
			relsN:=[ Vector(GF(p),b)*mat : b in relsN ];
			rels:=[];
			for rel in relsN do 				//hier falsche Bij. zwischen UF und F/Fneu???????????
				r:=Id(F);
				for j in [1..#Generators(F)] do
					r:=r*F.j^(Integers()!rel[j]);
				end for;
				Append(~rels,r);
			end for;
			//rels:=[F.1*F.2,F.3];
			Falt:=ncl<F|Fneu,rels>;
			G:=quo<F|Falt>;		//+ zus. Relation bei p=2???
			G;
			K:=L;
			_,_,E:=Factorization(Polynomial(K,s[2]):Extensions:=true);
			L:=E[1]`Extension;				//Faktor egal??
			Aut:=Automorphisms(K,CoefficientRing(f));
		else


			//Aut:=Automorphisms(K,CoefficientRing(f));  	//siehe unten!!
			normgrps:=[];
			for aut in Aut do 				//EZS würde auch reichen
				M:=CopyNormGroup(U,mU,N,aut);		//Kopieren in UF auch ausreichend!!!	
				Append(~normgrps,M);
			end for;
			N:=Intersec(normgrps);

			phi:=Inverse(mU)*mUF;                      
			Op1:=[];
			for i in [2..#Aut] do
				op:=[ Eltseq( phi(Aut[i](Inverse(phi)(UF.j))) ) : j in [1..n] ];
				Append(~Op1, GL(n,p)!Matrix(op));
			end for;
			Op1:=MatrixGroup<n,GF(p)|Op1>; //Operation auf der max el-ab Erw. aus der KKT (Autom)	
			//Frat,mFrat:=quo<Op1|FrattiniSubgroup(Op1)>;
			//Frat;
			//neuOp1:=MatrixGroup<n,GF(p)|[Inverse(mFrat)(Frat.k) : k in [1..#Generators(Frat)]]>;
			G,mG:=quo<F|Fneu>;
			//G;
			A,mA:=ncl<G|Generators(Falt)>;              	//NT entspricht max el-ab Erweiterung
			M,mM:=GModule(G,A,p);
			Op2:=Image(GModuleAction(M));
			//MM2:=GModule(Op2);
			//MM1:=GModule(neuOp1);
			//bool,B:=IsIsomorphic(MM1,MM2);
			//Op1;Op2;
			bool,B:=IsConjugate(GL(n,p),Op1,Op2);		//FEEEEEEEEHLER bei Q_8/D_4. Warum??
			if not bool then
				return bool;
			end if;
			//B:=B^-1;
			relsN:=[ mUF(U!N.i) : i in [1..#Generators(N)] ];
			relsN:=[ Eltseq(a) : a in relsN];
			relsN:=[ Vector(GF(p),a)*B : a in relsN ];
		
			rels:=[];
			for v in relsN do
				r:=Id(F);
				for i in [1..Dimension(M)] do 
					r:=r*Inverse(mG)(mA(Inverse(mM)(M.i)))^(Integers()!v[i]);   //Inverse(mG)(mA(Inverse(mM)(M.1))) : M -> G -> F
				end for;
				Append(~rels,r);
			end for;
			
			Fneu:=ncl<F|Fneu,rels>;
			G:=quo<F|Fneu>;
			IdentifyGroup(G);
			Falt:=Fneu;

			if i ne #s then						//im letzten Schritt keine explizite Konstr. des Körpers mehr nötig
				cfs:=ClassField(mU,N);     
				L:=CFComposite(cfs);
				
				AA:=[map<Integers(K)->Integers(L)|x:->aut(x)> : aut in Aut];				//Automorphisms() fehlerhaft
				cont:=[];										//besser selbst Autom. fortsetzen!!
				for i in [1..#AA] do
					AA[i]`LocalAutoGenerators:=ChangeUniverse(Aut[i]`LocalRingMap`LocalAutoGenerators, Integers(L));
					cont cat:=[ConstructFieldMap(L,m) : m in Continuations(AA[i],Integers(L))];
				end for;
				Aut:=cont;

				K:=L;
				_,_,E:=Factorization(Polynomial(K,s[i+1]):Extensions:=true);
				L:=E[1]`Extension;									//Faktor egal??
			end if;
		end if;
	end for;	
			
	return G;
end intrinsic;




intrinsic MyContinuations(K::FldPad, L::FldPad, aut::Map) -> SeqEnum
	{}
	a:=map<Integers(K)->Integers(L)|x:->aut(x)>;	
	cont:=[];
	a`LocalAutoGenerators:=ChangeUniverse(aut`LocalRingMap`LocalAutoGenerators, Integers(L));
	cont cat:=[ConstructFieldMap(L,m) : m in Continuations(a,Integers(L))];
	return cont;
end intrinsic;




intrinsic Test(f::RngUPolElt) -> SeqEnum//, SeqEnum
	{}
	mats:=[Matrix(GF(2),[[a1,a2,a3],[b1,b2,b3],[0,1,0]]) : a1,a2,a3,b1,b2,b3 in [0,1]];
	imats:=[ M : M in mats | IsInvertible(M) ];
	res:=[**];
	i:=0;
	for m in imats do
		i:=i+1;
		G:=FPGaloisGroup(f,m);
		if Type(G) eq GrpFP then
			Append(~res,<i,G,m>);
		end if;
	end for;
	return res;
end intrinsic;
		




//////////auch über Safarevic-Gruppe möglich://///////////////

intrinsic TSGaloisGroup(f::RngUPolElt) -> GrpFP             
	{Bestimmt die Galoisgruppe einer p^2-Erweiterung mit zwei Segmenten im RP. Zur Zeit nur p-Gruppen!!}
	
	//TO DO:
	//- Abfrage 2 Segmente, Test auf abelsch am Anfang? (KKT)
	//- zahmen Teil einbauen
	//- Option "Prec" bei UnitGroup. Dafür:
	//- Für die Operation auf dem el-ab NT immer die JNF nehmen? Normrelationen mit BWM umrechnen?

	Kf:=CoefficientRing(f);
	p:=Prime(Kf);
	s:=RamPolySubfields(f);
	K:=ext<Kf|Expand(s[2])>;
	//K:=ext<Kf|g>;
	A:=Remove(Automorphisms(K,Kf),1);

	U,iso:=UnitGroup(K);         //: Prec:=UnitPrecision(Integers(K)) liefert falsche Gruppe (s.u.)
	//U;
	UF,nhom:=quo<U|p*U>;
	//UF;
	n:=#Generators(UF);
	phi:=Inverse(iso)*nhom;

	_,_,E:=Factorization(Polynomial(K,f):Extensions:=true); //oder direkt aus RamPolyFactors?
	L:=E[1]`Extension;	

	N,Nmap:=NormGroup(L,iso);
	
	//im:={@ nhom(Nmap(N.i)) : i in [1..#Generators(N)] @};
	//for aut in A do 
	//	M,m:=CopyNormGroup(U,iso,N,aut);
	//	{@ nhom(m(M.i)) : i in [1..#Generators(M)] @};
	//	im:=im meet {@ nhom(m(M.i)) : i in [1..#Generators(M)] @};
	//end for;

	for aut in A do
		M:=CopyNormGroup(U,iso,N,aut);
		N:=N meet M;
	end for;
	//N;
	
	im:=[ nhom(U!N.i) : i in [1..#Generators(N)] ];	   //Erzeuger der Normgruppe des ZerfKp in UF
	//im;

	op:=[ Eltseq( phi(A[1](Inverse(phi)(UF.i))) ) : i in [1..n] ];  //Operation von Cp auf dem el.-ab. NT
	//op;

	F:=FreeGroup(n+1);
	rels:=[];
	for i in [1..n] do 
		Append(~rels, F.i^p=Id(F));
		wi:=Id(F);
		for j in [1..n] do
			wi:=wi*F.j^op[i][j];
		end for;
		Append(~rels, F.i^F.(n+1)=wi);
	end for;
	for k in [1..n] do
		for l in [k+1..n] do
			Append(~rels, F.k^F.l=F.k);
		end for;
	end for;
	Append(~rels, F.(n+1)^p=F.1);      	////auch bei p=2 korrekt? Begründung?? (Caputo/Safarevic nur für p>2)////
                                                ////nur korrekt, wenn F.1^F.(n+1)=F.1. (vgl. Pazderski)              ////
                                                ////ist z.b. bei UnitGroup mit Prec:=... nicht immer der Fall!       //// 
                         
	H:=quo<F|rels>;                        //maximale el.-ab. Erweiterung über Cp-Erw.

	rels2:=[];
	for s in [1..#im] do		      //zusätzliche Relationen aus der Normgruppe	
		us:=Id(H);
		for t in [1..n] do
			us:=us*H.t^Eltseq(im[s])[t];
		end for;
		Append(~rels2,us);
	end for;
	//rels2;
	return	H, quo<H|rels2>;
end intrinsic;	





intrinsic CopyPolynomial(f::RngUPolElt, aut::Map) -> RngUPolElt
	{}
	l:=[ aut(Coefficient(f,j)): j in [0..Degree(f)]];
	return Polynomial(CoefficientRing(f),l);
end intrinsic;



intrinsic CComposite(poly::SeqEnum) -> FldPad					//einfacher in RngLocA!!
	{}
	N:=ext<CoefficientRing(poly[1])|Expand(poly[1])>;   //müsste das Polynom aus dem TK-Turm sein, daher eisenstein!! nicht reduzibel?????
	for i in [2..#poly] do
		i;
		_,_,E:=Factorization(Polynomial(N,poly[i]):Extensions:=true);
		N:=E[1]`Extension;
	end for;
	return N;
end intrinsic;



intrinsic CGaloisGroup(f::RngUPolElt) ->  FldPad, SeqEnum
	{Nur einfaches Kopieren der Polynome ohne Klassenkörpertheorie!}
	
	Kf:=CoefficientRing(f);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;
	A:=Automorphisms(K,Kf);                
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	

	for i in [2..#s] do
		KX<x>:=PolynomialRing(K);
		1;
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		1;
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        //immer richtiger Faktor??
		gi;
		3;
		L:=Ki;
		for aut in Remove(A,1) do 				//ein EZS würde auch reichen!!
			hi:=CopyPolynomial(gi,aut);
			hi;
			hi:=Polynomial(L,hi);
			if not HasRoot(hi) then                               //Bedingung korrekt??
				4;
			        _,_,EE:=Factorization(hi:Extensions:=true);
		                L:=EE[1]`Extension;
				4;
			end if;
		end for;
		3;
		5;
		AA:=[map<Integers(K)->Integers(L)|x:->aut(x)> : aut in A];
		cont:=[];
		for i in [1..#AA] do
			AA[i]`LocalAutoGenerators:=ChangeUniverse(A[i]`LocalRingMap`LocalAutoGenerators, Integers(L));
			cont cat:=[ConstructFieldMap(L,m) : m in Continuations(AA[i],Integers(L))];
		end for;
		5;
		A:=cont;
		K:=L;	
	end for;
    
	return K,A;

end intrinsic;




intrinsic LocAGaloisGroup(f::RngUPolElt) ->  FldPad
	{Kopieren der Polynome, aber Gleichheitstest über NormGroup und Konstruktion in RngLocA.
	 Funktioniert nicht. "molten bug"!!}
	
	Kf:=CoefficientRing(f);
	s:=Reverse(RamPolySubfields(f));
	K:=ext<Kf|Expand(s[1])>;

	for i in [2..#s] do
		KX<x>:=PolynomialRing(K);
		1;
		_,_,E:=Factorization(Polynomial(K,s[i]):Extensions:=true);
		1;
		Ki:=E[1]`Extension;
		gi:=DefiningPolynomial(E[1]`Extension);        //immer richtiger Faktor??
		
		U,iso:=UnitGroup(K);
		N:=NormGroup(Ki,iso);

		2;
		A:=Automorphisms(K,Kf);                        //Autom`Berechnung bei Körpertürmen fehlerhaft!??
		if #A ne Degree(K,Kf) then
			error "Fehler in der Automorphismenberechnung!";
		end if;	
		2;
		3;
		L:=Ki;
		for aut in Remove(A,1) do 				//ein EZS würde auch reichen!!
			hi:=Expand(CopyPolynomial(gi,aut));
			Ti:=ext<K|hi>;
			4;
			M:=NormGroup(Ti,iso);
			4;
			if M ne N then                               //Bedingung korrekt??
				L:=LocalField(L,Polynomial(L,hi));
			end if;
		end for;
		3;
	        K:=RamifiedRepresentation(L);	
	end for;
    
	return K;

end intrinsic;





///////////////////////////old////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*


intrinsic RamPolyFactors(f::RngUPolElt) -> SeqEnum
	{Calculates the factorization of f corresponding to the ramification polygon.}

	L<alpha>:=TotallyRamifiedExtension(CoefficientRing(f),f);
	P<x>:=PolynomialRing(L);
	rho:=Evaluate(f,alpha*x+alpha) div (x*alpha^Degree(f));	
	v:=Vertices(NewtonPolygon(rho));
	factors:=[];

	for i in Reverse([2..#v]) do
		N:=Degree(rho);
		m:=Integers() ! v[i-1][1];
		E:=Integers() ! (v[i][1]-v[i-1][1]);
		H:=Integers() ! (v[i-1][2]-v[i][2]);
		e:=E div GCD(E,H);
		h:=H div GCD(E,H);
		
		LL<beta>:=ext<L|x^e-alpha>;		
		R<y>:=PolynomialRing(LL);
		srho:=Evaluate(rho,y*beta^h)/beta^(N*h);                       //"flatten" the last segment
		c:=[ Coefficient(rho,j)*beta^(h*(j-N)) : j in [m..N] ];
		srho2:=Polynomial(LL,c);
		srho1:=y^m;
		11;
		hl:=HenselLift(Polynomial(Integers(LL),srho),[srho1,srho2]);
		11;
		rho1:=Polynomial(L,Evaluate(hl[1],y/beta^h)*(beta^h)^m);
		rho2:=Polynomial(L,Evaluate(hl[2],y/beta^h)*(beta^h)^(N-m));    //should be over L (unique factor cor. to last segment)
		Append(~factors,rho2);
		
		rho:=rho1;
	end for;

	factors:=Reverse(factors);
	factors:=[ Evaluate(rhoi,(x-alpha)/alpha)*alpha^Degree(rhoi) : rhoi in factors ];    //back to f...
	factors[1]:=factors[1]*(x-alpha);

	return factors;

end intrinsic;








intrinsic RamPolyFactors(f::RngUPolElt) -> SeqEnum
	{Vorraussetzung: Factorization liefert die Faktoren in der durch das RP vorgegebenen Reihenfolge!}
	
	RP:=RamificationPolygon(f);
	
	L:=TotallyRamifiedExtension(CoefficientRing(f),f);
	P<x>:=PolynomialRing(L);
	fak:=Factorization(Polynomial(L,f));
	fak:=[a[1]:a in fak];
	//fak;
	
	list:=[]; Append(~list,x^0); prodlist:=list[1];
	j:=1; g:=x^0;

	for i in [1..#fak] do
		g:=g*fak[i];
		if Degree(g) eq RP[j+1][1]+1 then
			Append(~list, g div prodlist);
			j:=j+1;
			prodlist:=prodlist*list[j];
		end if;		
	end for;
	
	Remove(~list,1);
	return list;
	
end intrinsic;





//BEISPIELE://///////////////////////////////////////////////////////////////////////

> Attach("/home/reh/greve/magma/rampoly.m");
> Attach("/home/reh/greve/magma/subfields.m");
F:=FreeGroup(2);
F1:=ncl<F|F.1^3,F.2,(F.1,F.2)>;
> F2:=GStrich(F1,3);
G2,mG2:=quo<F|F2>;
> Order(G2);
A,mA:=ncl<G2|F.1^3,F.2,(F.1,F.2)>;
M,mM:=GModule(G2,A,3);






*/

