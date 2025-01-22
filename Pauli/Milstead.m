

function Swap(A,B)
	return B,A;
end function;


function Difference_XGCD(m,n)

//Input:  Two positive integers m,n
//Output: returns integer g and natural numbers x and y such that g=gcd(m,n) and g=xm-yn.

	g,x,y:=XGCD(m,n);
	state:=(x gt 0) and (y lt 0);
	if state then
		return x, -y;
	end if;
	while (state eq false) do
		x:=x+n;
		y:=y-m;
		state:=(x gt 0) and (y lt 0);
	end while;
	return x, -y;
end function;

function Segmental_Inertia(A)
	//A is already in ResidueClassField's PolynomialRing
	ff:=Factorization(A);
	degs:=[Degree(ff[i][1]): i in [1..#ff]];
	f:=LCM(degs);
	return f;
end function;

function UnityIndex(K,ei)
	//Finds [K(zeta_ei):K], i.e. smallest natural number f so that ei divides (-1 +q^f).
	V,n:=ResidueClassField(K);
	q:=#V;
	f:=1;
	state:=IsDivisibleBy(-1 + q^f, ei);
	while state eq false do
		f:=f+1;
		state:=IsDivisibleBy(-1 + q^f, ei);
	end while;
	return f;
end function;

function My_SubField(T,m)
	phi:=DefiningPolynomial(T);
	n:=Degree(phi);
	if false eq IsDivisibleBy(n,m) then
		error "Error: m does not divide the degree of T.";
	end if;
	C:=ConstantCoefficient(phi);
	U:=BaseRing(T);
	Ux<x>:=PolynomialRing(U);
	psi:=x^m + C;
	return ext<U|psi>;
end function;

function MyIsConjugate(H,K:S:=0)
	// W:=Galois_Starting_Group(phi);
	//MyIsConjugate(Alt(14),CyclicGroup(14):S:=W);
	if Degree(K) ne Degree(H) then
		error "Error:  Inputted Groups have different Degrees";
	end if;
	if Type(S) eq RngIntElt then
		n:=Degree(H);
		return IsConjugate(Sym(n),H,K);
	end if;
	return IsConjugate(S,H,K);
end function;

function ChangeGroupDegree(G,n)
//Input:  A permuation group G and a natural number n.
//Output:  A subgroup of Sn that is isomorphic to G.
	if Degree(G) eq n then
		return G;
	end if;
	_,Grp:=RegularRepresentation(G,sub<G|>);
	if Degree(Grp) eq n then
		return Grp;
	else
		T:=TransitiveGroups(n);
		S:=[t: t in T| (Order(t) eq Order(G))];
		state:=false;
		i:=1;
		if IsEmpty(S) then
			return "fail";
		end if;
		while state eq false do
			Candidate:=S[i];
			state:=IsIsomorphic(G,Candidate);
			i:=i+1;
		end while;
		return Candidate;
	end if;
end function;


function OptimizedUnramifiedExtensions(K,f)
// Input: K is an unramified extension of Qp or Zp, f is a natural number.
// Output: Returns a unramified extension L of K of degree f in absolute 
//         representation (that is as an extension over Zp or Qp) and a map from K to L.
	F:=InertiaDegree(K);
	Qp:=BaseRing(K);
	L:=UnramifiedExtension(Qp,f*F);
	rho:=MinimalPolynomial(K.1);
	Lx<x>:=PolynomialRing(L);
	rho:=Lx!rho;
	gamma_prime:=Roots(rho:Max:=1)[1][1];
	omega:=map<K->L | x:-> &+[(gamma_prime)^(i-1) *Coefficient(x,i): i in [1..F]]>;
	return L,omega;
end function;

function Make_Cyclotomic(K)
// Input:  K is an unramified extension of Qp or Zp.
// Output: The cyclotomic Unramified Extension L of Qp (Zp) 
//         so that [L:Qp]=[K:Qp] and a map tau from K to L.
	f:=InertiaDegree(K);
	Qp:=BaseRing(K);
	L:=CyclotomicUnramifiedExtension(Qp,f);
	rho:=MinimalPolynomial(K.1);
	Lx<x>:=PolynomialRing(L);
	rho:=Lx!rho;
	gamma_prime:=Roots(rho:Max:=1)[1][1];
	tau:=map<K->L | x:-> &+[(gamma_prime)^(i-1) *Coefficient(x,i): i in [1..f]]>;
	return L,tau;
end function;  


function IsThereaU(W)
	if true eq IsUnramified(W) then
		a:=false eq IsTotallyRamified(W);
	else
		while false eq IsUnramified(W) do
			W:=BaseRing(W);
		end while;
		a:=false eq IsTotallyRamified(W);
	end if;
	return a;
end function;


function RecoverR(T)
	if true eq IsThereaU(T) then
		U:=BaseRing(T);
		p:=Prime(T);
		h:=DefiningPolynomial(T);
		c:=ConstantCoefficient(h);
		w:=-c div p;
		V,n:=ResidueClassField(U);
		if Log(n(w)) eq 0 then
			r:=0;
		else
			r:=Log(n(w))/Log(V.1);
		end if;
	else
		p:=Prime(T);
		h:=DefiningPolynomial(T);
		c:=ConstantCoefficient(h);
		w:=-c div p;
		Zp:=PrimeRing(T);
		q:=Precision(Zp);
		V,n:=ResidueClassField(Zp);
		b:=(PrimitiveElement(V)@@n)^(p^q);
		if Log(n(w)) eq 0 then
			r:=0;
		else
			r:=Log(n(w))/Log(n(b));  //Denominator should be one.  Drop it?
		end if;
	end if;
	return r;
end function;


function ChangeT(t)
	if true eq IsThereaU(t) then
		U:=BaseRing(t);
		e:=RamificationIndex(t);
		r:=RecoverR(t);
		r:=Integers()!r;
		pi:=Prime(t);
		Ux<x>:=PolynomialRing(U);
		h:=x^e - pi*(U.1)^r;	
		T:=TotallyRamifiedExtension(U,h);
	else
		e:=RamificationIndex(t);
		Zp:=PrimeRing(t);
		p:=Prime(Zp);
		q:=Precision(Zp);
		r:=RecoverR(t);
		r:=Integers()!r;
		pi:=p;
		V,n:=ResidueClassField(Zp);
		b:=(PrimitiveElement(V)@@n)^(p^q);
		Zpx<x>:=PolynomialRing(Zp);
		h:=x^e - pi*(b)^r;
		T:=TotallyRamifiedExtension(Zp,h);
	end if;
	return T;
end function;


function EqualDegree_Case(T1,T2)
// Compositum of the two tamely ramified extensions T1 and T2
	e:=RamificationIndex(T1);
	gamma1:=ConstantCoefficient(DefiningPolynomial(T1));
	gamma2:=ConstantCoefficient(DefiningPolynomial(T2));
	U:=BaseRing(T1);
	V,n:=ResidueClassField(U);
	Vz<z>:=PolynomialRing(V);
	tau:=z^e-n(gamma2/gamma1);
	ff:=Factorization(tau);
	degs:=[Degree(ff[i][1]): i in [1..#ff]];
	f:=LCM(degs);
	return UnramifiedExtension(T1,f);
end function;


function CoPrimeDegree_Case(T1,T2)
// Compositum of the two tamely ramified extensions T1 and T2
	g:=DefiningPolynomial(T1);
	h:=DefiningPolynomial(T2);
	c:=ConstantCoefficient(h);
	d:=ConstantCoefficient(g);
	n:=Degree(g);
	m:=Degree(h);
	r,a,b:=Xgcd(m,n);
	Tx<x>:=PolynomialRing(T1);
	Cons:=(c/d)^b;
	ChangePrecision(~Cons,Precision(T1));
	Cons:=T1!Cons;
	return TotallyRamifiedExtension(T1,x^m - T1.1*Cons);
end function;


function Restructure_Tower(C,a)

	if a eq 1 then
		U:=C;
		if Degree(U) eq 1 then
			return BaseRing(U);
		end if;
		T1:=BaseRing(U);
		K:=BaseRing(T1);
		Zp:=PrimeRing(K);
		e:=RamificationIndex(T1);
		pi_K:=UniformizingElement(K);
		psi:=DefiningPolynomial(T1);
		gamma1:=-ConstantCoefficient(psi) div pi_K;
		if InertiaDegree(K,Zp) gt 1 then
			U_prime, omega:=OptimizedUnramifiedExtensions(K,InertiaDegree(U));			
			Uy<y>:=PolynomialRing(U_prime);			
			gamma1_prime:=omega(gamma1);
			pi_U_prime:=UniformizingElement(U_prime);
			phi:=y^e -gamma1_prime*pi_U_prime;
			phi:=Uy!phi;
			T:=TotallyRamifiedExtension(U_prime,phi);
			return T;
		else  //Inertia Degree eq 1.
			new_inertia:=InertiaDegree(U);
			U_prime:=UnramifiedExtension(Zp,new_inertia);
			Uy<y>:=PolynomialRing(U_prime);
			gamma_prime:=U_prime!gamma1;
			pi_U_prime:=UniformizingElement(U_prime);
			phi:=y^e - gamma_prime*pi_U_prime;
			phi:=Uy!phi;
			T:=TotallyRamifiedExtension(U_prime,phi);
			return T;
		end if;
	end if;  //omega in both cases is a map from K to the base ring of T.
	if a eq 2 then
		T2:=C;
		T1:=BaseRing(T2);
		K:=BaseRing(T1);
		deg:=RamificationIndex(T1)*RamificationIndex(T2);
		deg:=Integers()!deg;
		Cons:=Norm(ConstantCoefficient(DefiningPolynomial(T2)));
		Kx<x>:=PolynomialRing(K);
		T:=TotallyRamifiedExtension(K, x^deg + Cons);
	end if;

	if a eq 3 then
		V:=C;
		T1_prime:=BaseRing(V);
		S1_prime:=BaseRing(T1_prime);
		U_prime:=BaseRing(S1_prime);

		e1:=RamificationIndex(V);
		e2:=RamificationIndex(T1_prime);
		e3:=RamificationIndex(S1_prime);

		C0:=ConstantCoefficient(DefiningPolynomial(V));
		Cons:=Norm(Norm(C0));

		Uz<z>:=PolynomialRing(U_prime);

		T:=TotallyRamifiedExtension(U_prime, z^(e1*e2*e3) + Cons);
	end if;
	return T;
end function;


function Third_Case(T1,T2)
// Compositum of the tamely ramified extensions T1 and T2
	e1:=Degree(DefiningPolynomial(T1));
	e2:=Degree(DefiningPolynomial(T2));

	K:=BaseRing(T1);
	pi_K:=UniformizingElement(K);

	m:=GCD(e1,e2);
	q:=e1/m;
	r:=e2/m;

	q:=Integers()!q;
	r:=Integers()!r;

	if q eq 1 then
		T1,T2:=Swap(T1,T2);
		e1,e2:=Swap(e1,e2);
		q,r:=Swap(q,r);
	end if;

	gamma1:=-ConstantCoefficient(DefiningPolynomial(T1)) div pi_K;
	gamma2:=-ConstantCoefficient(DefiningPolynomial(T2)) div pi_K;

	Zp:=PrimeRing(T1);

	S1:=My_SubField(T1,m);
	S2:=My_SubField(T2,m);

	pi_S1:=UniformizingElement(S1);		
	pi_S2:=UniformizingElement(S2);

	U:=EqualDegree_Case(S1,S2);
	S1_prime:=Restructure_Tower(U,1);
	pi_S1_prime:=UniformizingElement(S1_prime);

	U_prime:=BaseRing(S1_prime);
	pi_U_prime:=UniformizingElement(U_prime);

	//To move things from K to U_prime

	f:=InertiaDegree(U_prime)/InertiaDegree(K);
	f:=Integers()!f;
	UU, omega:=OptimizedUnramifiedExtensions(K,f);	//omega maps K to UU;
	eps:=map<UU->U_prime|x:->&+[(U_prime.1)^(i-1)*Coefficient(x,i): i in [1..Degree(UU)]]>;  //compensate for magma not realizing U_prime=UU
	My_EmbeddingMap:=omega*eps;	//Mapping from: RngPad: K to RngPad: U_prime  
	//Find gamma1_prime for T1/S1
	RHS_1:=(-1)^(m+1)*gamma1*pi_K div Norm(pi_S1);
	RHS_1:=Integers()!RHS_1;
	if RHS_1 eq 1 then
		gamma1_prime:=1;
	elif (RHS_1 eq -1) and (IsOdd(m) eq true) then
		gamma1_prime:=-1;
	else
		VK,nk:=ResidueClassField(K);
		Vx<x>:=PolynomialRing(VK);
		tau_1:=x^m - nk(RHS_1);
		tau_1:=Vx!tau_1;
		state:=Roots(tau_1:Max:=1)[1][1];
		gamma1_prime:=(state @@ nk);
	end if;
	//Form T1_prime as extension of S1_prime
	Sz<z>:=PolynomialRing(S1_prime);
	T1_prime:=TotallyRamifiedExtension(S1_prime, z^q - gamma1_prime*pi_S1_prime);
	if r eq 1 then
		return T1_prime, 2;
	end if;
	//Find gamma2_prime for T2/S2
	RHS_2:=(-1)^(m+1)*gamma2*pi_K div Norm(pi_S2);
	RHS_2:=Integers()!RHS_2;
	if RHS_2 eq 1 then
		gamma2_prime:=1;
	elif (RHS_2 eq -1) and (IsOdd(m) eq true) then
		gamma2_prime:=-1;
	else
		VK,nk:=ResidueClassField(K);
		Vx<x>:=PolynomialRing(VK);
		tau_2:=x^m - nk(RHS_2);
		tau_2:=Vx!tau_2;

		state:=Roots(tau_2:Max:=1)[1][1];
		gamma2_prime:=(state @@ nk);
	end if;

	//Write pi_S2 in terms of pi_S1_prime.
	Cons:=My_EmbeddingMap(gamma2*pi_K);
	Cons:=S1_prime!Cons;
	temp:=z^m - Cons;

	pi_S:=Roots(temp:Max:=1)[1][1];

	//Form T2_prime as extension of S1_prime	
	gamma2_prime:=My_EmbeddingMap(gamma2_prime);
	gamma2_prime:=S1_prime!gamma2_prime;
	ChangePrecision(~pi_S, Precision(S1_prime));

	T2_prime:=TotallyRamifiedExtension(S1_prime, z^r -gamma2_prime*pi_S);

	T:=CoPrimeDegree_Case(T1_prime,T2_prime);
	a:=3;	
	return T,a;
end function;


function PairComposite(T1,T2)
// Copositum of tamely ramified extensions T1 and T2
	p1:=Prime(T1);
	p2:=Prime(T2);

	p1:=Integers()!p1;
	p2:=Integers()!p2;

	if p1 ne p2 then
		error "Extensions are not defined over the same p-adic ring";
	end if;

	K:=BaseRing(T1);
	Zp:=PrimeRing(T1);

	n:=Degree(DefiningPolynomial(T1));
	m:=Degree(DefiningPolynomial(T2));

	Comp:=0; //should never happen, just here to guarantee it runs.

	if m eq n then
		Comp:= EqualDegree_Case(T1,T2);
		a:=1;
		Comp:=Restructure_Tower(Comp,a);
		return Comp;
	end if;

	if GCD(m,n) eq 1 then
		Comp:=CoPrimeDegree_Case(T1,T2);
		a:=2;
		Comp:=Restructure_Tower(Comp,a);
		return Comp;
	end if;

	if (GCD(m,n) gt 1 and m ne n) then
		Comp,a:=Third_Case(T1,T2);
		Comp:=Restructure_Tower(Comp,a);
		return Comp;
	end if;
	return Comp;
end function;


function TameCompositum(T)
//Input a list of Tamely Ramified Extensions over a common local field K.
//Output Composite of all the extensions, taken two at a time.

	if #T eq 1 then
		return T[1];
	end if;

	if #T eq 2 then
		return PairComposite(T[1],T[2]);
	end if;


	Zp:=PrimeRing(T[1]);

	if #T gt 2 then
		C:=PairComposite(T[1],T[2]);
		count:=2;		//number of extensions used

		while count lt #T do
			Tame_ext:=T[count+1];
			Zp:=PrimeRing(Tame_ext);

			new_inertia:=InertiaDegree(Tame_ext,Zp);   //under Tame_ext;
			old_inertia:=InertiaDegree(C,Zp);         //under C, i.e. from previous loop
			
			if new_inertia lt old_inertia then
				e:=RamificationIndex(Tame_ext);

				U_prime:=BaseRing(C);
				K:=BaseRing(Tame_ext);

				pi_U_prime:=UniformizingElement(U_prime);
				pi_K:=UniformizingElement(K);

				psi:=DefiningPolynomial(Tame_ext);
				Uy<y>:=PolynomialRing(U_prime);

				gamma1:=-ConstantCoefficient(psi) div pi_K;

				// Move gamma1 to U_prime
				
				f:= old_inertia / new_inertia;
				f:=Integers()!f;
				UU, omega:=OptimizedUnramifiedExtensions(K,f);
				tau:=map<UU->U_prime|x:->&+[(U_prime.1)^(i-1)*Coefficient(x,i): i in [1..Degree(UU)]]>;

				gamma1_prime:=omega(gamma1);  //move to UU.

				gamma1_prime:=tau(gamma1_prime);  //move to U_prime.

				phi:=y^e -gamma1_prime*pi_U_prime;
				phi:=Uy!phi;

				TT:=TotallyRamifiedExtension(U_prime,phi);
			else
				TT:=Tame_ext;

			end if;

			C:=PairComposite(C,TT);
			count:=count+1;

		end while;
		return C;
	end if;
end function;


function TestPair(T1,T2)

	C:=Composite(T1,T2);
	e1:=RamificationIndex(C,PrimeRing(T1));
	f1:=InertiaDegree(C,PrimeRing(T1));

	P:=PairComposite(T1,T2);
	e2:=RamificationIndex(P,PrimeRing(T1));
	f2:=InertiaDegree(P,PrimeRing(T1));

	differ:=(e1-e2)^2 + (f1-f2)^2;

	return differ;


end function;



function Max_Tame_Subextension(phi)
vprint Galois,1:"GaloisGroupMilstead: Max_Tame_Subextension";
//An implementation of Greve-Pauli Theorem 9.1.
//That is a function that given an Eisenstein polynomial returns the maximal tamely ramified subfield of its splitting field.

//Input: an Eisenstein Polynomial of degree n=e0*p^w over Qp or an unramified extension of Qp

	if IsEisenstein(phi) eq false then
		error "Polynomial is not Eisenstein.";
	end if;

	K:=CoefficientRing(phi);		
	
	n:=Degree(phi);
	p:=Prime(PrimeRing(K));

	m:=Valuation(n,p);
	e0:=n div p^m;
	e0:=Integers()!e0;

	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n); //Ramification Polynomial
	rho:=Lx!rho;
	ChangePrecision(~rho,Precision(L));

	ramification_polygon := NewtonPolygon(rho);

	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1];

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	a:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices
	b:=[Integers()!vertices[i][2]: i in [1..#vertices]];	//list of y-coordinates of vertices
	
	if e0 eq 1 then
		l:=#slopes;
	else  //e0>1
		l:=#slopes-1;		//l+1 st segment is horizontal on x-axis
	end if;

	s:=[];

	for j in [1..(l+1)] do
		s[j]:=Valuation(a[j],p);
	end for;

	pi_L:=UniformizingElement(L);
	
	e:=[];
	h:=[];

	for i in [1..l] do
		e[i]:=Denominator(-slopes[i]);		//list of (negative) slope denominators
		h[i]:=Numerator(-slopes[i]);		//list of (negative) slope numerators
	end for;

	gamma:=[];
	
	A:=[];						//Residual Polynomials of segments
	f:=[];						//Segmental Inertia of segments

	RL,omega:=ResidueClassField(L);
	Ry<y>:=PolynomialRing(RL);			//Polynomial Ring in variable y over Residue Class Field of L
	
	RK,sigma:=ResidueClassField(K);
	Rs<r>:=PolynomialRing(RK);

	for i in [1..l] do
		di:=a[i+1]-a[i];
		fin:=Integers()!(di/e[i]);
		Ai:=Polynomial(Ry,[&+[omega(Coefficient(rho,j*e[i]+a[i]) div pi_L^(-j*h[i]+b[i]))*y^j: j in [0..fin]]]);
		Ai:=Ry!Ai;
		A[i]:=Ai;

		f[i]:=Segmental_Inertia(A[i]);		
	end for;

	Indicess:=[UnityIndex(K,e0*e[i]): i in [1..l]];
	F:=LCM(Flat([f,Indicess]));

	I:=UnramifiedExtension(K,F);    

	RI,tau:=ResidueClassField(I);
	Rw<w>:=PolynomialRing(RI);

	for i in [1..l] do

		Ai:=A[i];
		Ai:=Rs!Ai;			//Might be unnecessary, magma knows that RL=RK

		Ai:=Rw!Ai;
		gamma_i:=Roots(Ai:Max:=1);
		gamma_i:=gamma_i[1][1];		//note gamma_i lives in RI
		gamma_i:=(gamma_i @@ tau);	//find pre_image of gamma_i, this is what we want
		Append(~gamma,gamma_i);

	end for;

	Iz<z>:=PolynomialRing(I);
	phi_0:=ConstantCoefficient(phi);

	//phi_0:=mi(phi_0);		//Moves (maps) phi_0 to I

	Tame_extensions:=[]; 

	for i in [1..l] do
		
		_, bi ,di:=Xgcd(h[i],e[i]);
		vi:=e0 * p^(m-s[i]) + n +1;		

		Ti:=TotallyRamifiedExtension(I, z^(e0*e[i]) - ((-1)^vi)* phi_0*gamma[i]^(bi*n) );
		Append(~Tame_extensions,Ti);
	end for;

	Towe:=TameCompositum(Tame_extensions);

	return Towe;	

end function;



function TameGaloisGenerators(T)

	if false eq IsTamelyRamified(T) then
		error "Error: Extension is not Tamely Ramified.";
	end if;

	U:=BaseRing(T);
	Zp:=PrimeRing(T);
	p:=Prime(T);
	q:=p;
	pi:=p;
	
	f:=InertiaDegree(U);
	e:=RamificationIndex(T);

	if false eq IsDivisibleBy(((q^f)-1),e) then
		error "Error: Extension is not normal.";
	end if;
	
	r:=RecoverR(T);
	r:=Integers()!r;

	k:=(r*(q-1))/e;

	l:=(-1+q^f)/e;

	k:=Integers()!k;
	l:=Integers()!l;

	state1:=U.1^k;	
	state2:=U.1^l;	
	
	g1:=map<U->U|x:->x>;
	g2:=map<U->U|y:->y^q>;
	
	h1:=map<T->T|x:->x*state2>;
	h2:=map<T->T|y:->y*state1>;

	s:=map<T->T|x:->&+[(h1(T.1))^(j-1)*(&+[(g1(U.1))^(i-1)*Coefficient(Coefficient(x,j),i): i in [1..f]]): j in [1..e]]>;
	
	t:=map<T->T|y:->&+[(h2(T.1))^(j-1)*(&+[(g2(U.1))^(i-1)*Coefficient(Coefficient(y,j),i): i in [1..f]]): j in [1..e]]>;

	return s,t,r,k;

end function;


function Slope_of_Factor(phi)

	NP:=NewtonPolygon(phi);

	m:=Slopes(NP)[1];
	return m;
end function;



function CompareSlopes(x,y)

	m1:=Slope_of_Factor(x);
	m2:=Slope_of_Factor(y);

	return m1-m2;

end function;


function MakeEisenstein(phi)
//Input:  a single variable polynomial phi.
//Output:  if phi is Eisenstein, it is returned.  Otherwise, the extension (totally ramified) generated by phi is created and the DefiningPolynomial is returned.

	if IsEisenstein(phi) then
		return phi;
	else
		_,_,m:=Factorization(phi:Extensions:=true);
		return DefiningPolynomial(m[1]`Extension);
	end if;

end function;



function eltseq(a)
  L :=[];
  R := Parent(a);
  F, m := ResidueClassField(R);
  for v in [1..Precision(R)] do
    b := m(a);
    Append(~L,b);
    a := (a-b@@m) div UniformizingElement(R);
  end for;
  return L;
end function;


function seqelt(R,L)

  pi := UniformizingElement(R);
  F, m := ResidueClassField(R);
  return &+[L[i]@@m*pi^(i-1):  i in [1..#L]];

end function;


function IsSurjective_Finite(phi)
	//phi is a map or it is a polynomial that we are treating as a map

	if Type(phi) eq RngUPolElt then

		K:=CoefficientRing(phi);

		if Type(K) eq FldFin then
			K:=Set(K);
			n:=#K;

			MyImage:={Evaluate(phi,k): k in K};
			return K subset MyImage;
			
		else
			error "Error: Polynomial must be over a finite field.";
		end if;

	end if;
	
	if Type(phi) eq Map then

		K:=Domain(phi);
		
		if Type(K) eq FldFin then
			K:=Set(K);
			n:=#K;

			MyImage:={phi(k): k in K};
			return K subset MyImage;
				
		else
			error "Error:  Map's Domain must be a finite field";
		end if;

	end if;

	return 0; 

end function;


function Hasse_Herbrand(RamPolygon,lambda)

	rho:=Polynomial(RamPolygon);
	n:=Degree(rho);
	L:=CoefficientRing(rho);
	p:=Prime(PrimeRing(L));

	m:=Valuation(n,p);
	e0:=n div p^m;
	e0:=Integers()!e0;

	vertices:= LowerVertices(RamPolygon);     

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
	end if;

	if e0 gt 1 then

		Exclude(~vertices,<n,0>);//Throw out the point (n,0)

	end if;

	VTerms:={(1/n)*(vertices[i][2] + lambda*vertices[i][1]): i in [1..#vertices]};

	return Min(VTerms);	

end function;


function reduce_poly(f)

// f must be eisenstein
vprint Galois,1:"GaloisGroupMilstead: reduce_poly in",f;
  o := CoefficientRing(Parent(f));
  n := Degree(f);
  fseq := Eltseq(f);
  j := Min([n*Valuation(o!n)] cat [n*(Valuation(o!i)+Valuation(fseq[i+1])-1)+i:i in [1..n-1]]);
  prec := Ceiling((n+2*j+1)/n);  //Krasner bound
  oo := ChangePrecision(o,prec);
vprint Galois,1:"GaloisGroupMilstead: reduce_poly out",Polynomial(o,Polynomial(oo,f)); 
  return Polynomial(o,Polynomial(oo,f)), prec; 
end function;


function MongeReduction(psi,KrasnerBound)

	K:=CoefficientRing(psi);
	Kx<x>:=PolynomialRing(K);
	n:=Degree(psi);

	L:=ext<K|psi>;
	pi_L:=UniformizingElement(L);

	alpha:=L.1;
	Ly<y>:=PolynomialRing(L);

	rho:=Evaluate(psi, alpha*y + alpha) div (alpha^n); //Ramification Polynomial
	rho:=Ly!rho;
	ChangePrecision(~rho,Precision(L));
	
	RamPolygon:=NewtonPolygon(rho);

	vertices := LowerVertices(RamPolygon);	
	slopes := Slopes(RamPolygon)[1..#vertices-1];

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	New_psi:=x^n;

	for i in [0..(n-1)] do
		C:=Coefficient(psi,i);
		Cs:=eltseq(C);

		m:=1;

		while ( (n + n*Hasse_Herbrand(RamPolygon,m) - i) / n  ) lt KrasnerBound do

			j:=(n + n*Hasse_Herbrand(RamPolygon,m) - i) / n;
			
			//if 0 eq ( (i-n*Hasse_Herbrand(RamPolygon,m)) mod n) then

			if Denominator(j) eq 1 then

				j:=Integers()!j;

				if (-m) in slopes then

					l:=Position(slopes,-m);
					k:=vertices[l][1];

					k:=Integers()!k;

					RL,omega:=ResidueClassField(L);
					Rz<z>:=PolynomialRing(RL);			//Polynomial Ring in variable y over Residue Class Field of L
	
					RK,sigma:=ResidueClassField(K);
					Rs<r>:=PolynomialRing(RK);

					dl:=vertices[l+1][1]-vertices[l][1];
					fin:=Integers()!(dl/Denominator(m));
		
					A:=Polynomial(Rz,[&+[omega(Coefficient(rho,Integers()!(u*Denominator(m)+vertices[l][1])) div pi_L^(Integers()!(-u*Numerator(m)+vertices[l][2])))*z^u: u in [0..fin]]]);
					A:=Rz!A;
					A:=Rs!A;

					//Sm:=(r^k)*A;
					//Sm;

					rho_m:=Evaluate(rho,alpha^m *y);
					contalpha:=Min([Valuation(r): r in Eltseq(rho_m)]);

					rho_m:=rho_m div alpha^contalpha;

					Sm:=Rz!rho_m;

					
					if IsSurjective_Finite(Sm) then

						Cs[j+1]:=0;			

					end if;

				else
					Cs[j+1]:=0;
				end if;

			end if;

			m:=m+1;

		end while; 

		New_psi:=New_psi + x^i * seqelt(K,Cs);

	end for;

	return New_psi;
	
end function;




function FactorofDegree(phi,n)
//Input:  a polynomial phi and a natural number n.
//Output: a factor of phi that has degree n.  If such a factor doesn't exist, an error message is returned to the user.
	n;

	if IsIrreducible(phi) then		//might not need.  might slow it down
		
		if Degree(phi) eq n then
			return phi;
		else
			error "Error: Polynomial is irreducible and doesn't have a factor of that degree.";
		end if;
	end if;

	Fac:=Factorization(phi);
	state:=false;
	
	for i in [1..#Fac] do
		i;
		candidate:=Fac[i][1];
		D:=Degree(candidate);
		D;		

		if Degree(candidate) eq n then
			return candidate;
		end if;
	end for;

	error "Error:  Polynomial doesn't have a factor of that degree.";

end function;


function Oneslope_onefactor(R)

//Input:  list of factors of ramification polynomial, each of which has a newton polygon that consists of one segment.

//Output: list of factors of ramification polynomial so that there is one factor for each segment of the ramification polygon.

	fin:=#R;

	count:=1;

	if fin eq count then
		return R;
	end if;

	S:=[];		//will contain the slopes of the factors.

	N:=[];		//will contain the new list of factors to be outputted	

	for i in [1..#R] do
		S[i]:=Slope_of_Factor(R[i]);
	end for;

	
	while count le fin do

		prod:=1;		
		
		state:=true;

		while state eq true do
			current:=R[count];
			prod:=prod*current;

			if count lt fin then
				state:=S[count] eq S[count+1];
			else
				state:=false;
			end if;

			count:=count+1;
			
		end while;

		Append(~N, prod);

	end while;

	return N;

end function;				




function RamificationPolygonFactors(phi)
//Input:  An Eisenstein polynomial phi in K[x]

//Output;  Factors of phi corresponding to the segments of the ramification polygon of phi.

	K:=CoefficientRing(phi);		
	n:=Degree(phi);
	ChangePrecision(~phi,Precision(K));

	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n); //Ramification Polynomial
	rho:=Lx!rho;
	ChangePrecision(~rho,Precision(L));
	
	ramification_polygon := NewtonPolygon(rho);

	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1]; 


	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	rho:=rho div x;

	Facc:=Factorization(rho);
	R:=[Facc[i][1]: i in [1..#Facc]];

	Sort(~R,CompareSlopes);				//Reorder the factors by segment

	if #slopes ne #R then

		R:=Oneslope_onefactor(R);		//How I make the Ri's match up with the proper segments.
	end if;


	F:=[];						//Will  contain the factors of phi.


	F1:=(x-alpha)*Evaluate(R[1],((x-alpha)/alpha))*alpha^(Degree(R[1]));
	F[1]:=Lx!F1;

	for i in [2..#slopes] do
		Fi:=alpha^(Degree(R[i])) * Evaluate(R[i],((x-alpha)/alpha));
		F[i]:=Lx!Fi;
	end for;

	return F;

end function;



function RamificationPolygonSubfields(phi:St:=0)
//Input:  An Eisenstein polynomial phi in K[x].

//Output:  The set {L_{l+1}, ... , L_1} of subfields of K[x]/(phi).  L_j's aren't ordered.

	K:=CoefficientRing(phi);		
	n:=Degree(phi);

	phi:=reduce_poly(phi);
	ChangePrecision(~phi,Precision(K));
	
	Kx<x>:=PolynomialRing(K);

	L:=ext<K|phi>;
	alpha:=L.1;

	F:=RamificationPolygonFactors(phi);
	F_0:=[];   //Constant coefficients of F

	for i in [1..#F] do
		F_0[i]:=ConstantCoefficient(F[i]);
	end for;

	U:={};

	for i in [1..(#F-1)] do
		mu:= MinimalPolynomial(&*[F_0[j]: j in [1..i]]);	  //Compute the minimal polynomial mu in K[x] of Constant terms up to F_0[i]
		mu:=Kx!mu;


		if IsEisenstein(mu) then
			mu:=reduce_poly(mu);
			ChangePrecision(~mu,Precision(CoefficientRing(mu)));
			U:=U join {ext<K|mu>};

		else
			_,_,m:=Factorization(mu:Extensions:=true);
			U:=U join {m[1]`Extension};

		end if;

	end for;
	

	return U;

end function;



function RamificationPolygonTower(phi)

//Input:  An Eisenstein polynomial phi in K[x].
 
//Output:  The set {Li, ..., L1} such that the extension L=K[x]/(phi) is the tower of extensions L contains L1 contains ... contains Li contains K.

	phi:=MakeEisenstein(phi);
	K:=CoefficientRing(phi);		
	n:=Degree(phi);

	phi:=reduce_poly(phi);	
	ChangePrecision(~phi,Precision(K));

	Kx<x>:=PolynomialRing(K);

	L:=ext<K|phi>;
	alpha:=L.1;

	F:=RamificationPolygonFactors(phi);

	i:=#F -1;   

	if i eq 0 then   
		return [L];
	end if;

	H:=&*[F[j]: j in [1..i]]; 
	H_0:= ConstantCoefficient(H);
	//ChangePrecision(~H_0,Precision(K));

	mu:=MinimalPolynomial(H_0);
	mu:=Kx!mu;

	//Form new extension E

	if IsEisenstein(mu) then
		mu:=reduce_poly(mu);
		ChangePrecision(~mu,Precision(CoefficientRing(mu)));
		E:=ext<K|mu>;
	else

		_,_,m:=Factorization(mu:Extensions:=true);
		E:=m[1]`Extension;  //new extension
	end if;

	//reduce extension's defining polynomial

	psi:=DefiningPolynomial(E);
	psi:=reduce_poly(psi);

	U:=BaseRing(E);
	ChangePrecision(~psi,Precision(U));
	E:=ext<U|psi>;
	Ey<y>:=PolynomialRing(E);
	
	phi:= Ey!phi;

	H:=FactorofDegree(phi,Degree(H));	//Could pick wrong factor.  Not likely.  Check this in case of error down the road.	

	return [E] cat RamificationPolygonTower(H);

end function;


function TameGaloisGroup(T)

	n:=Degree(T, PrimeRing(T));
	
	if RamificationIndex(T,PrimeRing(T)) eq 1 then
		f:=n;

		return CyclicGroup(f);
	end if;


	if IsWildlyRamified(T) then
		error "Error:  Extension is not tame";
	end if;


	// Part 1: Get T in the right form

	if IsThereaU(T) eq false then
		T:=ChangeT(T);		//if I is degree 1, gets rid of it.
		

	else
		V:=BaseRing(T);
		f:=InertiaDegree(V);
		
		UU,mu:=Make_Cyclotomic(V);
		psi:=DefiningPolynomial(T);

		new_Coeffs:=[];

		for i in [0..Degree(psi)] do
			old_coeff:=Coefficient(psi,i);		
			Append(~new_Coeffs, mu(old_coeff));					
		end for;

		Towe:=TotallyRamifiedExtension(UU,Polynomial(UU,new_Coeffs));
		T:=ChangeT(Towe);
		
	end if;


	e:=RamificationIndex(T);
	f:=InertiaDegree(T,PrimeRing(T));

	r:=RecoverR(T);
	r:=Integers()!r;
	q:=#ResidueClassField(PrimeRing(T));

vprint Galois,1:"GaloisGroupMilstead: r", r;
vprint Galois,1:"GaloisGroupMilstead: f", f;

	

	//Part 2:  Check Normality.  Call bool.

	bool:=IsDivisibleBy(((q^f)-1),e) and IsDivisibleBy((r*(q-1)),e);

	

	//Part 3:  If bool is true.  

	if bool and (f eq 1) then
		return CyclicGroup(e);
	end if;


	if bool and (e eq 1) then

		return CyclicGroup(f);
	end if;


	if bool then

		G<x,y>:=Group<x,y|x^e, y^f=x^r, x*y=y*x^q>;  //x=s, y=t
		F:=f;

	end if;

	
	
	//Part 4:  If bool is false, form new extension that is normal

	if bool eq false then

		g:=GCD(q^f -1, r*(q-1));

		u:=1;
		RHS:=(e*(q^f - 1)/g);
		RHS:=Integers()!RHS;

		state:=IsDivisibleBy(((q^(f*u))-1),RHS);

		while state eq false do

			u:=u+1;
			state:=IsDivisibleBy(((q^(f*u))-1),RHS);
			u; //Maybe get rid of?

		end while;


		s:=r*((q^(f*u))-1) div (q^f -1);
		s:=Integers()!s;

		G<x,y>:=Group<x,y|x^e, y^(f*u)=x^s, x*y=y*x^q>;
		
	end if;


	PP,mp:=PermutationGroup(G);	

	return ChangeGroupDegree(PP,n);


end function;




function  One_Segment(phi)
vprint Galois,1:"GaloisGroupMilstead: One Segment";
//Input: phi Eisenstein in OK[x] of degree p^m such that the ramification polygon of phi has one segment.

//Output:  Gal(phi) as a subgroup of AGL_m(F_p).

	
	if IsEisenstein(phi) eq false then
		error "Polynomial is not Eisenstein.";
	end if;

	K:=CoefficientRing(phi);
	n:=Degree(phi);
		
	if Type(K) eq FldPad then

		error "Coefficient Ring must be a ring, not a field.";
	end if;
	
	Fp:=ResidueClassField(PrimeRing(K));
	p:=#Fp;
	m:=Valuation(n,p);

	Fq:=ResidueClassField(K);
	q:=#Fq;

	L:=ext<K|phi>;
	alpha:=L.1;
	Lx<x>:=PolynomialRing(L);
	pi_L:=UniformizingElement(L);

	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n); //Ramification Polynomial
	rho:=Lx!rho;
	ChangePrecision(~rho,Precision(L));

	Segmentslope:=Slope_of_Factor(rho div x);

	e:=Denominator(-Segmentslope);
	h:=Numerator(-Segmentslope);

	RL,omega:=ResidueClassField(L);
	Ry<y>:=PolynomialRing(RL);

	Fx:=PolynomialRing(Fq);

	di:=n-1;
	fin:=Integers()!(di/e);

	bi:=h*(n-1) div e;
	A:=Polynomial(Ry,[&+[omega(Coefficient(rho,j*e+1) div pi_L^(-j*h+bi))*y^j: j in [0..fin]]]);
	A:=Ry!A;
	A:=Fx!A;		//Residual Polynomial

	f1:=Segmental_Inertia(A);
	f:=LCM(f1,UnityIndex(K,e));

	b, b_tilde:=Difference_XGCD(h,e);
	a, a_tilde:=Difference_XGCD(e,n);

	Fqf:=ResidueClassField(UnramifiedExtension(K,f));
	Fqf_unitgroup, unitmap:= UnitGroup(Fqf);
	zeta:=unitmap(Fqf_unitgroup.1);

	Fz<z>:=PolynomialRing(Fqf);
	A:=Fz!A;

	//Now to find the roots of A and r.

	RR:=Roots(A);
	u:=[r[1]: r in RR];
	u1:=u[1];

	if Log(u1) eq 0 then
		r_prime:=0;
	else
		r_prime:=Integers()!(Log(u1^b)/Log(zeta));
	end if;

	
	if r_prime in  Seqset([0..(e-1)]) then
		r:=r_prime;

	else
		r:=(r_prime mod e);

	end if;	


	//Finding M

	Fqf_plus,ADDMap:=AdditiveGroup(Fqf); //ADDMap maps Fqf_plus to Fqf
	M:=sub<Fqf_plus|Identity(Fqf_plus)>;
	i:=1;

	state:=Order(M) eq p^m;


	while state eq false do

		Cons:=u[i]/(zeta^(r*h));
		Cons:=Fqf!Cons;
		psi:=z^e -Cons;

		RRi:=Roots(psi);
		ui:=[r[1]: r in RRi];
		gens:=[(a*s)@@ADDMap: s in ui];
		M:=sub<Fqf_plus|gens cat Setseq(Generators(M))>;
		i:=i+1;

		state:=Order(M) eq p^m;

	end while;
	
	MGens:=Generators(M);
	V:=VectorSpace(Fp,#MGens);
	B:=Basis(sub<V|[V!Eltseq(b): b in MGens]>);

	k:=(r*(q-1))/e;
	l:=(-1+q^f)/e;

	k:=Integers()!k;
	l:=Integers()!l;

	Expss:=[];	

	for j in [1..Dimension(V)] do

		current:=M.j;
		current:=ADDMap(current);

		if Log(current) eq 0 then
			Expss[j]:=0;
		else
			Expss[j]:=Integers()!(Log(current)/Log(zeta));

		end if;				

	end for;


	S_rows:=[];
	T_rows:=[];


	for s in [1..Dimension(V)] do		//Finding the rows of S using the automorphism s defined on powers of zeta

		j:=Expss[s];

		temp:=zeta^(l*h +j);
		temp:=temp@@ADDMap;
		temp:=M!temp;
		temp:=V!Eltseq(temp);

		temp:=ElementToSequence(temp);

		S_rows[s]:=temp;

	end for;


	for t in [1..Dimension(V)] do		//Find the rows of T using the automorphism t defined on powers of zeta

		j:=Expss[t];

		temp:=zeta^(h*k +q*j);
		temp:=temp@@ADDMap;
		temp:=M!temp;
		temp:=V!Eltseq(temp);

		temp:=ElementToSequence(temp);

		T_rows[t]:=temp;  		

	end for;


	S:=Matrix(GF(p), S_rows);
	T:=Matrix(GF(p), T_rows);

	GGLL:=GL(m,Fp);

	S:=GGLL!S;
	T:=GGLL!T;

	Matgenn:=sub<GGLL|S,T>;

	G:=AffineGroup(Matgenn);
	
	return G;
		
	//return S,T;

	
end function;


function My_IsomorphismClasses(lis)

	//A:=AssociativeArray(Sym(Degree(lis[1])));   //Associative array with index universe S_n.
	A:=AssociativeArray();
	GroupOrders:={Order(g): g in lis};

	GroupsbySize:={};   //Each element is a set of all the groups of a particular order
	lis:=Set(lis);

	for n in GroupOrders do

		G:={};
		for H in lis do
			if Order(H) eq n then
				G:=G join {H};
				lis:=lis diff {H};
			end if;
		end for;

		GroupsbySize:=GroupsbySize join {G};
	end for;


	for m in GroupsbySize do

		collection:=m;

		while #collection gt 0 do
			Candidate:=Random(collection);
			IClass:={Candidate};
			collection:=collection diff IClass;

			for g in collection do
				if IsIsomorphic(g,Candidate) then
					IClass:=IClass join {g};
					collection:=collection diff {g};
				end if;
			end for;
	
			A[Candidate]:=IClass;			

		end while;

	end for;

	return A;

end function;


function TransitiveIDClasses(lis)

	A:=AssociativeArray();
	GroupOrders:={Order(g): g in lis};

	GroupsbySize:={};   //Each element is a set of all the groups of a particular order
	lis:=Set(lis);

	for n in GroupOrders do

		G:={};
		for H in lis do
			if Order(H) eq n then
				G:=G join {H};
				lis:=lis diff {H};
			end if;
		end for;

		GroupsbySize:=GroupsbySize join {G};
	end for;

	for m in GroupsbySize do

		collection:=m;

		while #collection gt 0 do
			Candidate:=Random(collection);
			IClass:={Candidate};
			collection:=collection diff IClass;

			for g in collection do
				if TransitiveGroupIdentification(g) eq TransitiveGroupIdentification(Candidate) then
					IClass:=IClass join {g};
					collection:=collection diff {g};
				end if;
			end for;
	
			A[Candidate]:=IClass;			

		end while;

	end for;

	return A;

end function;


function CandidateClasses(lis)

	G:=lis[1];
	n:=Degree(G);

	if n le TransitiveGroupDatabaseLimit() then
		return TransitiveIDClasses(lis);
	else
		return My_IsomorphismClasses(lis);
	end if;

	return 0;

end function;

// load "AbsoluteResolvents.magma";


function My_Degree(A)
	if (Type(A) eq RngUPolElt) then
		r:=Degree(A);
	else
		r:=0;
	end if;

	if A eq 0 then
		r:=0;
	end if;

	return r;

end function;



function IsitSquareFree(f)

	if Type(CoefficientRing(f)) eq FldFin then
		f:=f;
	else

		Z:=Integers();
		Zx<x>:=PolynomialRing(Z);
		f:=Zx!f;
	end if;

	answer:=My_Degree(GCD(f,Derivative(f))) eq 0;

	return answer;

end function;


function Tschirnhausen(P)

	if (Type(P) ne RngUPolElt) then
		error "Error: P is not univariate.";
	end if;

	Z:=Integers();
	Zw<w>:=PolynomialRing(Z);
	P:=Zw!P; //Comment out if statement if you want to work over ring other than the integers.
	
	R:=CoefficientRing(P);
	T:=[];
	n:=My_Degree(P);
	
	S:=Eltseq(P);

	Rxy<x,y>:=PolynomialRing(R,2);
	
	Py:=Evaluate(P,y);	
	Py:=Rxy!Py;

	Q:=x^2;
	state:=UnivariatePolynomial(GCD(Q,Derivative(Q,x)));  //Note:  unlike most of my code, state is not a boolean value.

	while ( state notin R) do		//My_Degree(GCD(f,Derivative(f))) gt 0;
		for y in [1..n] do
			a:=Random(0,50);
			Append(~T,a);
		end for;
		
		Ay:=Polynomial(Rxy,T);
		Ay:=Evaluate(Ay,y);
	
		v:=Rxy!(x-Ay);
		
		Q:=Resultant(Py,v,y);

		state:=UnivariatePolynomial(GCD(Q,Derivative(Q,x)));

		if My_Degree(state) eq 0 then
			state:=ConstantCoefficient(state);
		end if;

		T:=[];

	end while;

		Rx<x>:=PolynomialRing(R);	
		Q:=UnivariatePolynomial(Q);
		Q:=Rx!Q;
		
	return Q;

end function;


function sz(f,g)
	//Input: polynomials f,g
	
	//Computes the resultant Res_y(f(y),g(x-y)).


	Z:=Integers();
	Zw<w>:=PolynomialRing(Z);
	f:=Zw!f;
	g:=Zw!g;

	R:=CoefficientRing(f);
	n:=Degree(f);
	Rxy<x,y>:=PolynomialRing(R,2);


	First:=Rxy!Evaluate(f,y);
	Second:=Rxy!Evaluate(g,(x-y));

	Res:=Resultant(First,Second,y);
	Res:=UnivariatePolynomial(Res);

	return Res;

end function;



function mz(d,f)
	//Input:  polynomial f, integer d
	//Output:  the polynomial whose roots are equal to the roots of f multiplied by the integer d.

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);
	f:=Zx!f;

	n:=Degree(f);
	
	if d eq 0 then

		return x^n;
	end if;

	//Need to evaluate f at y/d so temporarily change polynomial ring so that I'm over the rationals.

	Q:=Rationals();
	Qy<y>:=PolynomialRing(Q);
	f:=Qy!f;

	goal:=(d^n)*Evaluate(f,y/d);
	goal:=Zx!goal;

	return goal;

end function;



function pr(k,u)

	//Input: Positive integer k, monic polynomial u such that u(x)=r(x)^k for some unknown monic r(x)

	//Output:  r(x)

	if k eq 1 then
		return u;
	end if;

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);
	u:=Zx!u;

	t:=u div GCD(u,Derivative(u));
	r:=t;
	s:=u;

	while Degree(r) lt Degree(u)/k do
		s:=s/(t^k);
		s:=Zx!s;		

		t:=GCD(s,t);
		r:=t*r;

	end while;

	return r;

end function;


function dp(f)

	//computes the familiar resolvent of degree n*(n-1)/2, corresponding to the invariant T = x_1 + x_2 that is stabilized by S_2 x S_{n-2}. 

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);
	f:=Zx!f;

	First:=sz(f,f);
	Second:=mz(2,f);

	First:=Zx!First;
	Second:=Zx!Second;

	return pr(2,Zx!(First/Second));

end function;



function lr2(f,a,b)

	//lr2(f,1,2) computes the resolvent of degree n*(n-1), corresponding to T = x_1 + 2*x_2 that is stabilized by S_1 x S_1 x S_{n-2}.
	// It should not matter what the last two arguments are, as long as they are non equal integers.  

	if a eq b then
		error "Error:  Last two arguments should not be equal.";
	end if;

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);

 	R:= sz(mz(a,f),mz(b,f))/mz(a+b,f);

	return (Zx!R);

end function;




function tp(f) 

	//computes the resolvent of degree n*(n-1)*(n-2)/6, corresponding to T = x_1 + x_2 + x_3 that is stabilized by S_3 x S_{n-3}.

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);


	R:= pr(3,sz(dp(f),f)/lr2(f,1,2));	

	return (Zx!R);

end function;



function lr112(f)

	// computes the resolvent of degree n*(n-1)*(n-2)/2, corresponding to T = x_1 + x_2 + 2*x_3 that is stabilized by S_2 x S_1 x S_{n-3}.


 	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);

	R:= sz(dp(f),mz(2,f))/lr2(f,1,3);

	return (Zx!R);

end function;


function qp(f) 

	//computes the resolvent of degree n*(n-1)*(n-2)*(n-3)/24, corresponding to T = x_1 + x_2 + x_3 + x_4 that is stabilized by S_4 x S_{n-4}.

	Z:=Integers();
	Zx<x>:=PolynomialRing(Z);

	R:=pr(4,sz(tp(f),f)/lr112(f));

	return (Zx!R);

end function;



function tally(seq)

	set := Seqset(seq);
	multiset:=SequenceToMultiset(seq);
	return Sort([<s,Multiplicity(multiset,s)>: s in set]);

end function;


/////Below functions also used for relative resolvents

function NeedResolvent(W,H,Candidates)

	if #Candidates eq 1 then
		return false;
	end if;

	a,b,c:=CosetAction(W,H);

	//Figure out orbitlengths for first candidate.  This will be the basis of comparison.

	G1:=Candidates[1];

	phi_G1:=a(G1);

	Orbb1:=Orbits(phi_G1);
	Orb_lengths1:=[#r: r in Orbb1];
	
	Orb_lengths1:=Set(tally(Orb_lengths1));

	for i in [2..#Candidates] do

		G:=Candidates[i];
		phi_G:=a(G);

		Orbb:=Orbits(phi_G);
		Orb_lengths:=[#r: r in Orbb];
	
		Orb_lengths:=Set(tally(Orb_lengths));
	
		if Orb_lengths ne Orb_lengths1 then
			return true;
		end if;

	end for;

	return false;

end function;



function dpFact(phi)

	K:=CoefficientRing(phi);
	p:=Prime(K);

	psi:=dp(phi);

	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=dp(phi);
		state:=IsitSquareFree(psi);

	end while;

	
	FF:=pFactorDegrees(p,psi);

	return Set(tally(FF));

end function;



function lrFact(phi)
	//Resolvent factorization for the lr2 resolvent described above.

	K:=CoefficientRing(phi);
	p:=Prime(K);

	psi:=lr2(phi,1,2);

	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=lr2(phi,1,2);
		state:=IsitSquareFree(psi);

	end while;

	
	FF:=pFactorDegrees(p,psi);

	return Set(tally(FF));

end function;



function tpFact(phi)

	K:=CoefficientRing(phi);
	p:=Prime(K);

	psi:=tp(phi);

	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=tp(phi);
		state:=IsitSquareFree(psi);

	end while;
	
	FF:=pFactorDegrees(p,psi);

	return Set(tally(FF));

end function;



function  S3_OrbitLengthTest(G,phi,ResolventFactorList)

	//G is a candidate
	//tp(f) computes the resolvent


	n:=Degree(phi);
	Sn:=Sym(n);

	H:=DirectProduct(Sym(3),Sym(n-3));   //S_3 x S_{n-3}
	m:=Order(Sn)/Order(H);
	m:=Integers()!m;

	a,b,c:=CosetAction(Sn,H);

	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	return Orb_lengths eq ResolventFactorList;


end function;


function lr_OrbitLengthTest(G,phi,ResolventFactorList)

	//G is a candidate
	
	n:=Degree(phi);
	Sn:=Sym(n);

	H:=DirectProduct(Sym(1),DirectProduct(Sym(1),Sym(n-2)));   //S_1 x S_1 x S_{n-2}
	m:=Order(Sn)/Order(H);
	m:=Integers()!m;

	a,b,c:=CosetAction(Sn,H);

	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	return Orb_lengths eq ResolventFactorList;

end function;



function LRFact(phi)

	//Resolvent factorization for the lr112 resolvent described above.

	K:=CoefficientRing(phi);
	p:=Prime(K);

	psi:=lr112(phi);
	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=lr112(phi);
		state:=IsitSquareFree(psi);

	end while;

	
	FF:=pFactorDegrees(p,psi);

	return Set(tally(FF));

end function;


function LR_OrbitLengthTest(G,phi,ResolventFactorList)

	//G is a candidate
	
	n:=Degree(phi);
	Sn:=Sym(n);
	
	H:=DirectProduct(Sym(2),DirectProduct(Sym(1),Sym(n-3)));   //S_2 x S_1 x S_{n-3}
	a,b,c:=CosetAction(Sn,H);

	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	return Orb_lengths eq ResolventFactorList;

end function;


function qpFact(phi)

	//Resolvent factorization for qp resolvent mentioned above

	K:=CoefficientRing(phi);
	p:=Prime(K);

	psi:=qp(phi);

	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=qp(phi);
		state:=IsitSquareFree(psi);

	end while;

	
	FF:=pFactorDegrees(p,psi);

	return Set(tally(FF));

end function;


function qp_OrbitLengthTest(G,phi,ResolventFactorList)

	//G is a candidate
	
	n:=Degree(phi);
	Sn:=Sym(n);
	
	H:=DirectProduct(Sym(4),Sym(n-4));  //S_4 x S_{n-4}
	a,b,c:=CosetAction(Sn,H);

	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	return Orb_lengths eq ResolventFactorList;

end function;


// load "Brown.magma";


function  IspGroup(G,p)

    return Order(G) eq p^(Valuation(Order(G),p));

end function;



function QuickIsSubgroup(H,G) 

//Smaller Group first

//Is H conjugate to a subgroup of G?  Doing a partial check rather than computing subgroups of G since that can be time consuming.

//Instead we check the order and exponent of the two groups as well as whether or not they are abelian and whether or not they are cyclic.  

	B:=[true];

	if IsDivisibleBy(Order(G),Order(H)) then

		Append(~B, true);
	else

		return false;
	end if;

	
	if IsDivisibleBy(Exponent(G),Exponent(H)) then

		Append(~B, true);
	else

		return false;
	end if;

	if IsAbelian(G) then

		if IsAbelian(H) then
			Append(~B, true);
		else
			return false;
		end if;
	end if;


	if IsCyclic(G) then

		if IsCyclic(H) then
			Append(~B, true);
		else
			return false;
		end if;
	end if;


	return Set(B) eq {true};

end function;


function GlobalGalFilter(Candidates,phi)

	n:=Degree(phi);

	Z:=Integers();
	Zw<w>:=PolynomialRing(Z);

	Zphi:=Zw!phi;

	W:=GaloisGroup(Zphi);

	if IsIsomorphic(Sym(n), W) then
		return Candidates;
	end if;

	if IsIsomorphic(Alt(n), W) then
		return Candidates;
	end if;


	RemainingCandidates:=[G: G in Candidates|QuickIsSubgroup(G,W)];

	return RemainingCandidates;


end function;


function Only_One(lis)

//Confirms that, up to isomorphism, there is only one candidate group left.

	G:=lis[1];
	i:=1;

	while i lt #lis do
		H:=lis[i+1];
		bool:=IsIsomorphic(G,H);

		if bool eq false then
			return false;
		end if;

		i:=i+1;
	end while;
		
	return true;

end function;


function halve(a)
	
	a:=a/2;
	a:=Round(a);
	return a;

end function;


function SequenceToTuple(lis)

	if Type(lis) ne SeqEnum then

		error "Error: input should be a sequence.";
	end if;

	len:=#lis;

	A:=<>;  //output

	count:=1;

	while count le len do

		Append(~A,lis[count]);
		count:= count +1;

	end while;

	return A;


end function;



function TupleToSequence(t)

	if Type(t) ne Tup then
		error "Error: input should be a tuple.";
	end if;

	len:=#t;

	A:=[];

	for i in [1..len] do
 
		A[i]:=t[i];
	end for;

	return A;

end function;



function Number_Of_Segments(phi)

	K:=CoefficientRing(phi);		
	Zp:=BaseRing(K);

	n:=Degree(phi);
	p:=Prime(PrimeRing(K));

	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n);
	rho:=Lx!rho;

	ramification_polygon := NewtonPolygon(rho);
	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1]; 

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	return #slopes;

end function;



function StemField(phi)

	Qp:=CoefficientRing(phi);

	if IsInertial(phi) then
		K:=ext<Qp|phi>;
	elif IsEisenstein(phi) then
		K:=ext<Qp|phi>;
	else
		_,_,m:=Factorization(phi:Extensions:=true);
		K:=m[1]`Extension;

		psi:=DefiningPolynomial(K);
		psi:=reduce_poly(psi);

		U:=BaseRing(K);
		ChangePrecision(~psi,Precision(U));

		KK:=ext<U|psi>;
		return KK;
	end if;

	return K;

end function;



function GroupParity(G)
 
   if G subset AlternatingGroup(Degree(G)) then
     return 1;
   else 
     return -1;
   end if;

end function;


function SegmentOrderSieve(candidates,order,p)

	Latest_candidates:=[];
	fin:=#candidates[1];

	for j in [1..#candidates] do
		G:=candidates[j][fin];
		LL:=LowIndexSubgroups(G,<order,order>);

		if #LL gt 0 then
			for L in LL do
				if (Exponent(quo<G|L>) eq p) and IsElementaryAbelian(quo<G|L>) then
					Append(~Latest_candidates,Append(candidates[j],L));
				end if;
			end for;
		end if;
	end for;

	return Latest_candidates;


end function;


function PolynomialParity(phi)

	discc:=Discriminant(phi);

	if IsSquare(discc) then
		return 1;
	else
		return -1;
	end if;

end function;


function LowDegreeResolventTest(G,p)

	if p in {2,3,5} then
		return true;
	end if;

	if IsSolvable(G) eq false then
		error "Error:  Group G is not solvable";
	end if;

	Q:=PrimesInInterval(3,p-1);	//Possible primes

	for q in Q do
		if IsDivisibleBy(p-1,q) then

			H:=LowIndexSubgroups(G,<q,q>);
				for h in H do
					// _,g,_:=CosetAction(G,h);
					g:=CosetImage(G,h);

					state:=IsCyclic(g);
					
					if state eq false then
						return false;
					end if;

				end for;

		end if;

	end for;

	return true;

end function;


function Possible_Inertia(G,e0,p,f)

//Note:  remove condition c lt Degree(G) when you make this faster.

	choices:=[];

	remain:=Order(G) div e0;    //Remove |G0:G1| from picture
	remain:=Integers()!remain;

	highestpow:=Valuation(remain,p);

	lowestinertia:= remain div (p^(highestpow)); 
	lowestinertia:= Integers()!lowestinertia;

	for i in [0..highestpow] do
		temp:=lowestinertia * p^i ;
		Append(~choices,temp);
	end for;

	final_choices:=[c: c in choices | IsDivisibleBy(c, f) and c lt Degree(G)];

	return final_choices;


end function;


function Ramification_Segments(stem)


	if IsThereaU(stem) eq true then   //mixed case.  Stem field is totally ramified extension over unramified extension

		psi:=DefiningPolynomial(stem); //Eisenstein
		U:=BaseRing(stem);
		zp:=PrimeRing(stem);

		Qp:=pAdicField(Prime(zp),Precision(zp));
		UU := UnramifiedExtension(Qp,Polynomial(Qp,DefiningPolynomial(U)));  
		Coeff:=Coefficients(psi);
		n:=Degree(psi);

		
		tau:=map<U->UU|x:->&+[(UU.1)^(i-1)*Coefficient(x,i): i in [1..Degree(U)]]>;
		Ux<x>:=PolynomialRing(UU);

		phi:=Polynomial(Ux,[&+[tau(Coeff[i])*x^(i-1): i in [1..(n+1)]]]);
		phi:=Ux!phi;



	else		//stem field is just totally ramified extension generated by eisenstein poly.

 
		psi:=DefiningPolynomial(stem);  //Eisenstein
		zp:=BaseRing(stem);

		Qp:=pAdicField(Prime(zp),Precision(zp));
		Qx<x>:=PolynomialRing(Qp);
		phi:=Qx!psi;

	end if;

	RR:=RamificationPolygonTower(phi);
	R:=[];			//what's to be returned.

	for i in [1..#RR] do
		R[i]:=RingOfIntegers(RR[i]);
	end for;

	//Reduce defining polynomial of top extension.
	
	RL:=R[#R];

	V:=BaseRing(RL);
	tau:=DefiningPolynomial(RL);
	tau:=reduce_poly(tau);

	ChangePrecision(~tau,Precision(V));
	R[#R]:=ext<V|tau>;

	return R;

end function;



function MultiSolvable(n,p)

    t := TransitiveGroups(n);
    sol := [a : a in t|IsSolvable(a)];
    // solvable transitive groups of degree n

    s:=[b: b in sol|LowDegreeResolventTest(b,p)];


    l := [LowIndexSubgroups(a,<1,n>):a in s];
    // subgroups of index 1 to n.  Maybe do 2 to n.  
 

    candidates := []; // possible Galois groups
    for i in [1..#l] do
         vprint Galois,1:"GaloisGroupMilstead: \n#G",#s[i];

          for g0 in l[i] do
                if IsNormal(s[i],g0) and IsCyclic(quo<s[i]|g0>) then
                      g0index := Index(s[i],g0);
                          if true then
                           vprint Galois,1:"GaloisGroupMilstead:   [G:G0]",g0index;
                            G0:= g0;
                            lg0 := LowIndexSubgroups(g0,<1,p^g0index-1>);
                           
                        for h in lg0 do
                                  vprint Galois,1:"GaloisGroupMilstead: [G0:G1]",Index(g0,h);
                                      if IsCyclic(quo<g0|h>) and (p^g0index-1) mod Index(g0,h) eq 0 then
                                            if IsNormal(s[i],h) and IspGroup(h,p) then   //Last one not included in Brown paper nor awtrey one.  Not needed since subgroup of solvable group is solvable.
                                              vprint Galois,1:"GaloisGroupMilstead:     [G0:G1]",Index(g0,h);
                                              G1:= h;
                                              Append(~candidates,<s[i],G0,G1>);
                                              vprint Galois,1:"GaloisGroupMilstead:     #candidates",#candidates;
                                        end if;
                                end if;
                        end for;
                    end if;
            end if;        
        end for;
    end for;

    return candidates;


end function;


function StemFieldStartingGroup(stem)

	prod:=Sym(1);	//initial factor

	RR:=Ramification_Segments(stem);
	L:=RR[#RR];	//tower of subfields of eisenstein part of stem field. 
	
	WreathFactors:=[]; 

	//when doing loop on L...compare L to bottom of tower.

	while L ne PrimeRing(L) do
				
		if IsWildlyRamified(L) then
			
			f:=DefiningPolynomial(L);
			temp:=One_Segment(f);
			L:=BaseRing(L);

		else
			temp:=TameGaloisGroup(L);
			L:=PrimeRing(L);

		end if;
			
		Append(~WreathFactors,temp);
				
	end while;

	WreathFactors:=Reverse(WreathFactors);

	for i in [1..#WreathFactors] do
		Factt:=WreathFactors[i];
		prod:=WreathProduct(Factt,prod);
	end for;

	return prod;

end function;


function StartingTransitiveGroups(stem)

	G:=StemFieldStartingGroup(stem);
	Can:=Subgroups(G:IsTransitive:=true);

	C:=[r`subgroup: r in Can];
	C:=CandidateClasses(C);

	return Keys(C);

end function;



function Galois_Starting_Group(phi)

	stem:=StemField(phi);
	
	return StemFieldStartingGroup(stem);

end function;



function StartingGroupTest(phi,Candidate_number)

	if IsInertial(phi) then

		return 	Candidate_number eq 1;	

	end if;

	stem:=StemField(phi);
	T:=StartingTransitiveGroups(stem);

	Candidates:=[TransitiveGroupIdentification(t): t in T];

	bool:=Candidate_number in Candidates;

	return bool;

end function;


function CorrectE(phi)


	//Recovering e0
	K:=CoefficientRing(phi);
	n:=Degree(phi);
	p:=Prime(PrimeRing(K));
	w:=Valuation(n,p);
	e0:=n div p^w;
	e0:=Integers()!e0;


	//Create Ramfication Polygon

	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n);
	rho:=Lx!rho;

	ramification_polygon := NewtonPolygon(rho);
	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1]; 
	
	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;
	
	//Find rest of e_i's

	if e0 eq 1 then
		l:=#slopes;
	else  //e0>1
		l:=#slopes-1;		//l+1 st segment is horizontal on x-axis
	end if;

	m:=slopes;			//list of slopes

	e:=[];

	for i in [1..#slopes] do
		e[i]:=Denominator(-m[i]);		//list of (negative) slope denominators
	end for;

	e:=e[1..l];

	return (e0*LCM(e)), vertices, slopes;


end function;



function Max_Tame_Test(Candidate, TameGrp, phi)

	if IsEisenstein(phi) eq false then
		error "Polynomial is not Eisenstein.";
	end if;


	if Degree(Candidate) ne Degree(phi) then
		error "You didn't use regular representation";
	end if;

	PP:=TameGrp;
	Zp:=PrimeRing(CoefficientRing(phi));
	p:=Prime(Zp);

	if Order(TameGrp) eq 1 then
		return IspGroup(Candidate,p);
	end if;	

	ord:=Order(PP);

	Pos:=LowIndexSubgroups(Candidate,<ord,ord>);

	if IsEmpty(Pos) then 
		return false, [0];
	end if;

	B:=[];  //holds the possible subgroups S if any.

	for S in Pos do 
		bool:=IsIsomorphic(PP, quo<Candidate|S>) and IsNormal(Candidate,S) and IspGroup(S,p);

		if bool then
			Append(~B,S);
		end if;
	end for;

	if #B gt 0 then
		return true, B;
	else
		return false, [0];
	end if;

	return 0;

end function;



function GrieveSieve(phi)

	n:=Degree(phi);

	if IsInertial(phi) then
		return {TransitiveGroupIdentification(CyclicGroup(n))};
	end if;

	Qp:=CoefficientRing(phi);
	p:=Prime(Qp);
	purity:=1;  //if stays 1 then the polynomial is either eisenstein or inertial


	if Valuation(n,p) eq 0 then

		TameGrp:=TameGaloisGroup(StemField(phi));
		return {TransitiveGroupIdentification(TameGrp)};
	end if;


	if n eq p then
		phi:=MakeEisenstein(phi);
		Galgrp:=One_Segment(phi);
		return {TransitiveGroupIdentification(Galgrp)};
	end if;


	if IsEisenstein(phi) then
		K:=ext<Qp|phi>;
			if (Number_Of_Segments(phi) eq 1) and (n eq (p^(Valuation(n,p)))) then
				Galgrp:=One_Segment(phi);
				return {TransitiveGroupIdentification(Galgrp)};
			end if;
	 
	else   //Later Fix:  it is possible for 
		_,_,ext:=Factorization(phi:Extensions:=true);
		K:=ext[1]`Extension;   //new_phi for the Max Tame?
		purity:=0;
	
	end if;


	if Valuation(Degree(DefiningPolynomial(K)),p) eq 0 then
		TameGrp:=TameGaloisGroup(K);
		return {TransitiveGroupIdentification(TameGrp)};
	end if;
		

	Autoo:=AutomorphismGroup(K,Qp);

	if Order(Autoo) eq n then
		return {TransitiveGroupIdentification(Autoo)};  //We are in the normal case here.
		
	end if;

	t := StartingTransitiveGroups(K);
	//t := TransitiveGroups(n);
	
	sol:=[k: k in t|(Order(k) gt n)];  //Remove groups of order n.
	m:=[b: b in sol|LowDegreeResolventTest(b,p)];
	//m:=sol;
	
	PParity:=PolynomialParity(phi);

	S:=[];

	for i in [1..#m] do
		G:=m[i];
		Cen:=Centralizer(Sym(n),G);
	
		if (PParity eq GroupParity(G)) and IsIsomorphic(Autoo,Cen) then
			Append(~S,G);
		end if;

	end for;

	if Only_One(S) then
		return {TransitiveGroupIdentification(S[1])};
	end if;
	
		
	if purity eq 0 then
		phi:=DefiningPolynomial(K);
		n:=Degree(phi);
	end if;

	_,vertices,slopess:=CorrectE(phi);
	T:=Max_Tame_Subextension(phi);
	e0:=RamificationIndex(T, PrimeRing(T));
	f:=InertiaDegree(T,PrimeRing(T));

	l:=[];

	for j in [1..#S]  do 	//loop over S add things to l one at a time
		G:=S[j];
		Possiblef:=Possible_Inertia(G,e0,p,f);
	
		l[j]:=Flat([[LowIndexSubgroups(G,<Possiblef[i],Possiblef[i]>): i in [1..#Possiblef]]]);
	end for;
	

	candidates := []; // possible Galois groups
		
	for i in [1..#l] do
 		"\n#G",#S[i];

	          for g0 in l[i] do
               		if IsNormal(S[i],g0) and IsCyclic(quo<S[i]|g0>) then
               			g0index := Index(S[i],g0);
                       			if true then
                       				"  [G:G0]",g0index;
                       				G0:= g0;
                       				lg0 := LowIndexSubgroups(g0,<e0,e0>);
                       				for h in lg0 do
                               				vprint Galois,1:"GaloisGroupMilstead: [G0:G1]",Index(g0,h);
                               				if IsCyclic(quo<g0|h>) and (p^g0index-1) mod Index(g0,h) eq 0 and (p^g0index-1) ge e0 then
                                       				if IsNormal(S[i],h) and IspGroup(h,p) then   //Last one not included in Brown paper nor awtrey one.  Not needed since subgroup of solvable group is solvable.
                                       					"    [G0:G1]",Index(g0,h);
                                       					G1:= h;
                                       					Append(~candidates,<S[i],G0,G1>);
                                       					"    #candidates",#candidates;
                                       				end if;
                               				end if;
                       				end for;
               				end if;
        		end if;        
        	end for;
    	end for;


	New_candidates:=candidates;
	

	MainGroups:=Setseq({New_candidates[i][1]: i in [1..#New_candidates]});  //kill redundancy in G.

	if Only_One(MainGroups) then
		return {TransitiveGroupIdentification(r): r in MainGroups};
	end if;


	//Now need to incorporate the elementary abelian extensions from x-coordinates of Ramification Polygon.

	state:=false;

	while state eq false do
	
		if Reverse(vertices)[1][1] gt p^(Valuation(n,p)) then 	//have a horizontal segment 
			vertices:=vertices[1..(#vertices-1)];  		// remove last vertex
		end if;

		state:=(Reverse(vertices)[1][1] eq p^(Valuation(n,p)));
	end while;


	count:=#vertices;

	segmentorders:=[]; //will hold the degrees of the extensions.

	if count eq 1 then   //no non-horizontal segments, no help here.  For now just return answers so far.  Will need to be fixed later.
		return {TransitiveGroupIdentification(r): r in MainGroups};

	else 	while count gt 1 do
			current:=vertices[count][1]/vertices[count-1][1];
			current:=Integers()!current;
			Append(~segmentorders,current);
			count:=count-1;
		end while;
	end if;


	for i in [1..#segmentorders] do

		New_candidates:=SegmentOrderSieve(New_candidates, segmentorders[i],p);

		if IsEmpty(New_candidates) then
			return {TransitiveGroupIdentification(r): r in MainGroups};
		end if;

		ans:=Setseq({New_candidates[i][1]: i in [1..#New_candidates]});
		ret:={TransitiveGroupIdentification(r): r in ans};

		if #ret eq 1 then

			return ret;
		end if;
	end for;

	MainGroups:=Setseq({New_candidates[i][1]: i in [1..#New_candidates]});  

	if Only_One(MainGroups) then
		return {TransitiveGroupIdentification(r): r in MainGroups};
	end if;

	//If Eisenstein use Max_Tame_Test(Candidate, TameGrp, phi)
	
	if purity eq 1 then
		TameGrp:=TameGaloisGroup(T);
		"burp";
		UptoTame:=[m: m in MainGroups|Max_Tame_Test(m,TameGrp, phi)];
		MainGroups:=Set(UptoTame);
	end if;			
	 
	return {TransitiveGroupIdentification(r): r in MainGroups};	

end function;


function RecoverSlopes(phi)

	if IsInertial(phi) then
		return 0;
	end if;

	if IsEisenstein(phi) eq false then
		_,_,m:=Factorization(phi:Extensions:=true);
		K:=m[1]`Extension;
		phi:=DefiningPolynomial(K);
	end if;
		

	//Recovering e0
	K:=CoefficientRing(phi);
	n:=Degree(phi);
	p:=Prime(PrimeRing(K));
	w:=Valuation(n,p);
	e0:=n div p^w;
	e0:=Integers()!e0;


	//Create Ramfication Polygon

	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n);
	rho:=Lx!rho;

	ramification_polygon := NewtonPolygon(rho);
	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1]; 
	
	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;
	
	return slopes,vertices;

end function;


function S2Resolvent(f)

	//Computes the Resultant (f(y), f(x+y))/x^n.  Resultant done in the variable y.  Last, halve the exponents.

	Z:=Integers();
	Zw<w>:=PolynomialRing(Z);
	f:=Zw!f;

	R:=CoefficientRing(f);
	n:=Degree(f);
	Rxy<x,y>:=PolynomialRing(R,2);

	f1:=Evaluate(f,y);
	f1:=Rxy!f1;

	f2:=Evaluate(f,x+y);
	f2:=Rxy!f2;

	Res:=Resultant(f1,f2,y);
	Res:=UnivariatePolynomial(Res);

	Res:= Res div UnivariatePolynomial(x^n);  //Fine to here

	deg:=Degree(Res);
	Coeff:= Coefficients(Res);

	Rx<x>:=PolynomialRing(R);
	temp:=&+[Coeff[i]*(x)^(halve(i-1)): i in [1..(deg+1)]];
	Res:=Rx!temp;

	return Res;

end function;



function  S2_OrbitLengthTest(G,phi,ResolventFactorList)

	//G is a candidate

	n:=Degree(phi);
	Sn:=Sym(n);

	H:=DirectProduct(Sym(2),Sym(n-2));
	m:=Order(Sn)/Order(H);
	m:=Integers()!m;

	a,b,c:=CosetAction(Sn,H);

	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	return Orb_lengths eq ResolventFactorList;


end function;


function PossibleGalois(phi)

	Candidates:=GrieveSieve(phi);
	
	if #Candidates eq 1 then
		return Candidates;
	end if;

	Candidates:={TransitiveGroup(Degree(phi),c): c in Candidates};

	K:=CoefficientRing(phi);	

	psi:=S2Resolvent(phi);

	state:=IsitSquareFree(psi);

	while (state eq false) do
		phi:=Tschirnhausen(phi);
		psi:=S2Resolvent(phi);
		state:=IsitSquareFree(psi);

	end while;

	p:=Prime(K);
	FF:=pFactorDegrees(p,psi);
	S2ResolventFactorList:=Set(tally(FF));

	New_Candidates:=[G: G in Candidates |S2_OrbitLengthTest(G,phi,S2ResolventFactorList)];

	WhatRemains:={TransitiveGroupIdentification(r): r in New_Candidates};

	return WhatRemains; 

end function;



//Stuff for confirming sinclair.

function My_RamificationPolygon(phi)

	if IsEisenstein(phi) eq false then
		error "Polynomial is not Eisenstein.";
	end if;

	K:=CoefficientRing(phi);		
	
	n:=Degree(phi);
	L:=ext<K|phi>;
	alpha:=L.1;

	Lx<x>:=PolynomialRing(L);
	rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n); //Ramification Polynomial
	rho:=Lx!rho;
	ChangePrecision(~rho,Precision(L));

	ramification_polygon := NewtonPolygon(rho);

	return ramification_polygon,rho;  

end function;


function Recover_Vertices_Slopes(ramification_polygon)

	vertices := LowerVertices(ramification_polygon);	
	slopes := Slopes(ramification_polygon)[1..#vertices-1];

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	return vertices, slopes;

end function;


function RecoverPolyDegree(lis)
	//Uses the ramification polygon's vertices to recover the degree of the polynomial.

	vertices:=lis;
	a:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices

	return Integers()!Max(a);

end function;



function RecoverDiscriminant(lis)
	//recover discriminant's valuation from first vertex of polygon given as list of vertices

	if lis[1][1] eq 0 then
		lis:=lis[2..#lis];		//get rid of point with x-coordinate=0
	end if;

	if lis[1][1] ne 1 then
		error "Error: polygon doesn't have a vertex with x-coordinate 1.";
	end if;

	j:=lis[1][2];
	n:=RecoverPolyDegree(lis);


	return n+j -1;

end function;


// load "AMP-Resolvent.magma";


function fac(f,p)

	Fp:=FiniteField(p);
	Fpx<x>:=PolynomialRing(Fp);
	ff:=Factorization(Fpx!f);
	
	return [<Degree(g[1]),g[2]>: g in ff];

end function;



function ChoosePrime(f)

	//p:=Degree(f);
	//np:=Max(20,p);  //number of primes to find.  Will be decreased.

	p:=1;
	np:=20 + Degree(f); //number of primes to find.  Will be decreased.

	lp:=[];

	d:=Discriminant(f)*LeadingCoefficient(f);

	//while (#lp lt 1) and (np gt 0) do

	while #lp lt 1 do

		while (np gt 0) do
			p:=NextPrime(p);
		
			if d mod p ne 0 then

				np:=np-1;
				h:=PolynomialRing(FiniteField(p))!f;

				ff:=Factorization(h);
				degs:=[Degree(ff[i][1]): i in [1..#ff]];

				Append(~lp, <p,Sort(degs),Lcm(degs)>);
			end if;

		end while;

		np:=10;

	end while;

	Third:=[];  //list of Lcm(degs) = inertia.

	for j in [1..#lp] do
		Third[j]:=lp[j][3];
	end for;

	m,pos:=Minimum(Third);

	q:=lp[pos][1];

	return q,m;

end function;



function PartialTschirnhausen(P,q)
        Z:=Integers();
        Zx<x>:=PolynomialRing(Z);
       
	n:=My_Degree(P);
        T:=[Random(0,q):i in [1..n-2]];
	A :=  x+x^2*Polynomial(Z,T);
vprint Galois,1:"GaloisGroupMilstead: PartialTschirnhausen",A;        
        return A;

end function;


function MyComplexRootss(f,Prec)

	f:=PolynomialRing(Integers())!f;

	C:= ComplexField(Prec);
	Cx<x>:=PolynomialRing(C);

	f:=Cx!f;

	RRoots:=Roots(f);
	alpha:=[];

	for i in [1..#RRoots] do
		alpha[i]:=RRoots[i][1];
	end for;

	return alpha;

end function;


function FujiwaraBound(f)

  //Returns the Fujiwara bound for the absolute value of the complex roots of f

	c := Eltseq(f); _n := #c;

	max_comp := Maximum([Abs(c[_n-i]/c[_n])^(1/i) : i in [1.._n-2]] cat [0]);
	max_comp := 2*Maximum(max_comp, Abs(c[1]/2/c[_n])^(1/(_n-1)));
   
	return Ceiling(max_comp); //Fujiwara bound

end function;



function MY_LiftRoots(GloGenerator,alpha)

//double precision of roots.

	precis:=Precision(alpha[1]);
	C:=Parent(alpha[1]);
	k:=2*precis;

	ChangePrecision(~C,k);
	Cx<x>:=PolynomialRing(C);
	//psi:=Cx!GloGenerator;
	psi:=GloGenerator;
	Z:=Integers();
	beta:=[];

	for j in [1..#alpha] do
		beta[j]:=HenselLift(Cx!psi,C!alpha[j],k);
		//beta[j]:=HenselLift(psi,C!alpha[j],k);
	end for;

	return beta;

end function;


function lift_roots(Psi,alpha,k)


	//Compare precision of alpha to k.  If it's already greater than or equal to k don't lift.

	Currentprec:=Precision(Parent(alpha[1]));

	if Currentprec ge k then
		return alpha;
	end if;


	//Otherwise, lift and maintain root order

	C:=Parent(alpha[1]);
	q:=Prime(C);   
        ChangePrecision(~C,k);
 	Cx<x>:=PolynomialRing(C);
	psi:=Cx!Psi;
	beta:=[];

	for j in [1..#alpha] do
		beta[j]:=HenselLift(Cx!psi,C!alpha[j],k);   
	end for;

        return beta;

end function;


function pAdicRoot_Approx(f,Precis)

	f:=PolynomialRing(Integers())!f;
	
	p,m:=ChoosePrime(f);
	Prec:=10 + Valuation(Discriminant(f),p) + Precis;

	Zp:=pAdicRing(p,Prec);
	
	if m gt 1 then
		U:=UnramifiedExtension(Zp,m);
	else
		U:=Zp;
	end if;

	Ux:=PolynomialRing(U);
	f:=Ux!f;

	RRoots:=Roots(f);
	alpha:=[];

	for i in [1..#RRoots] do
		alpha[i]:=RRoots[i][1];
	end for;

	return alpha;

end function;



function VerifyInvariant(F,G,H)
	//Check that F is a G-relative H-invariant.
	//Currently uses the IsInvariant command which is true to high probability.  May replace this with commented bit if not convinced.

	if (H subset G) eq false then

		error "Error: H is not a subgroup of G.";
	end if;

	if Index(G,H) eq 1 then
		error "Error:  H equals G.";
	end if;

	n:=Degree(G);

	G:=Set(G);
	H:=Set(H);
	D:=G diff H;
	
	//Hcheck:={IsInvariant(F,sigma):sigma in H};  //Correct with high probability.
	//Gcheck:={IsInvariant(F,sigma):sigma in D};  //Check difference in G and H


	//Hcheck:={(F^sigma eq F):sigma in H};
	//Gcheck:={(F^sigma eq F):sigma in D}; //Not wise.  ONly recognizes equality when sigma is identity.


	Z:=Integers();
	Zx<[X]>:=PolynomialRing(Z,n);
	xvector:=[];

	for i in [1..n] do

		xvector[i]:=X[i];

	end for;
	
	Hcheck:={(Evaluate(F^sigma,xvector) eq Evaluate(F,xvector)):sigma in H};  //Not a good idea for high degree but irrefutable.
	Gcheck:={(Evaluate(F^sigma,xvector) eq Evaluate(F,xvector)):sigma in D}; 

	boolH:=Hcheck eq {true};

	if boolH  then

		boolG:=Gcheck eq {false};

		if boolG then

			return true;

		end if;

	end if;

	return false;

end function;



function Resolvent(G,H,pRootss,M,GloGenerator)

//M is the complex bound on the pRoots.
	Z:=Integers();

	if (H subset G) eq false then

		error "Error: H is not a subgroup of G.";
	end if;

	if Index(G,H) eq 1 then

		error "Error:  H equals G.";
	end if;

	//Establish F, a G-relative H-invariant	

	F:=RelativeInvariant(G,H);

	//if VerifyInvariant(F,G,H) then
		vprint Galois,1:"GaloisGroupMilstead: Invariant PASSES";
	//else
		vprint Galois,1:"GaloisGroupMilstead: Invariant FAILS";
	//end if;

	CReps:=RightTransversal(G,H);
	alpha:=pRootss;
        
	N:=Bound(GaloisGroupInvariant(G,H),Abs(M));  //to be modified later

	C:=Parent(alpha[1]);
	q:=Prime(C);   
	k:=400*Ceiling(Log(q,2*N+1) +1);
	k;
	
        beta := lift_roots(GloGenerator,alpha,k);

	C:=Parent(beta[1]);
	ChangePrecision(~C,k);
	Cx<x>:=PolynomialRing(C);

	modulus := Prime(C)^Precision(C);

	// Now find coset reps for G/H

	ResolventFactors:=[];

	for sigma in CReps do
		Cons:=Evaluate(F^(sigma),beta);
		Append(~ResolventFactors,x - C!Cons);
	end for;

	Res:=&*ResolventFactors;
	Res:=Cx!Res;

	coeff:=[];

	for j in [0..(Degree(Res))] do
		a := Z!Coefficient(Res,j);
 		if a lt (modulus div 2) then
 			coeff[j+1]:=a;
 		else 
                	coeff[j+1]:=a-modulus;
		end if;
 	end for;

	RF:=Polynomial(Z,coeff);

	if IsitSquareFree(RF) then
		return RF;
	end if;

	state:=false;

	while (state eq false) do
		A:=PartialTschirnhausen(GloGenerator,q); 
		"Tschirn";

		N:=Bound(GaloisGroupInvariant(G,H),Abs(Evaluate(A,M)));
		k:=400*Ceiling(Log(q,2*N+1) +1);
		k;
		
		beta := lift_roots(GloGenerator,beta,k);

		L:=Parent(beta[1]);
		ChangePrecision(~L,k);
		Ly<y>:=PolynomialRing(L);

		A:=Ly!A;
		gamma:=[];

		for k in [1..#beta] do
			gamma[k]:=Evaluate(A,beta[k]);
		end for;

		ResolventFactors:=[];
		modulus := Prime(L)^Precision(L);

		for sigma in CReps do
 			Cons:=Evaluate(F^(sigma),gamma);
 			Append(~ResolventFactors,(y - L!Cons));  //Should it be x?  Shouldn't i be using Ly or is it the same as Cx?
		end for;

		Res:=&*ResolventFactors;
		Res:=Cx!Res;

		coeff:=[];

 		for j in [0..(Degree(Res))] do
			a := Z!Coefficient(Res,j);
 			if a lt (modulus div 2) then
 		                coeff[j+1]:=a;
 		        else 
                                coeff[j+1]:=a-modulus;
                        end if;
 		end for;

 		RF:=Polynomial(Z,coeff);
		state:=IsitSquareFree(RF);
        	
	end while;

	return RF;

end function;


function GGlobalRepresentation(phi,L0)

	phi:=reduce_poly(phi);

	K:=CoefficientRing(phi);
	n:=Degree(phi);

	KK:=PolynomialRing(K);
	phi:=KK!phi;

	ChangePrecision(~phi,Precision(K));

	Q:=Rationals();
	Z:=Integers();
	Zw<w>:=PolynomialRing(Z);

	R:=Z;
	Rxy<x,y>:=PolynomialRing(R,2);
	
	L0x:=PolynomialRing(L0);
	L1:=BaseRing(L0);
	g1:=DefiningPolynomial(L1,PrimeRing(L1));

	g1:=Zw!g1;
	g1:= Evaluate(g1,y);
	g1:=Rxy!g1;		//Maybe comment out.

	psi:=DefiningPolynomial(L0);
	psi:=reduce_poly(psi);
	ChangePrecision(~psi,Precision(L1));

	f:=&+[x^i*(&+[y^(j-1)*Integers()!Coefficient(Coefficient(psi,i),j): j in [1..Degree(g1)]]): i in [0..Degree(psi)]];
	f:=Rxy!f;

	GloGenerator:=Resultant(f,g1,y);
	GloGenerator:=UnivariatePolynomial(GloGenerator);

	bool:=IsIrreducible(KK!GloGenerator) and IsIsomorphic(phi,KK!GloGenerator);

	while bool eq false do

vprint Galois,1:"GaloisGroupMilstead: Bad GloGenerator";
		t:=PartialTschirnhausen(phi,50);	//Let t(y) be a random polynomial
		t:=Evaluate(t,y);
		t:=Rxy!t;

		f:=Evaluate(f, [x-t,y]);      //Let the new f be f(x-t(y),y) 
		f:=Rxy!f;

		GloGenerator:=Resultant(f,g1,y);
		GloGenerator:=UnivariatePolynomial(GloGenerator);	//Then compute the resultant of this new f and the original g.
		bool:=IsIrreducible(KK!GloGenerator) and IsIsomorphic(phi,KK!GloGenerator);  

	end while;

	H:=ext<Q|GloGenerator>;
	Hr<r>:=PolynomialRing(H);

	g1:=UnivariatePolynomial(g1);	
	RR:=Roots((Hr!g1):Max:=1)[1][1];
	
	h1:=PolynomialRing(Q)!Eltseq(RR);  //Embedding polynomial

	m1:=Degree(g1);
	Block_Length:=(Integers()!(n/m1));

	First:=ext<Q|g1>;
	Fz<z>:=PolynomialRing(First);
	
	ff:=Factorization(Fz!GloGenerator);

	tau:=ff[1][1];
	//Wholee:=ext<First|tau>;

	Gal1:=GaloisGroup(tau);
	Gal2:=GaloisGroup(g1);

	Wprime:=WreathProduct(Gal1,Gal2);
	ba, im, ker:=BlocksAction(Wprime,{1..Degree(Gal1)});
	WreathBlocksprime:=Orbits(ker);

	WreathBlocksprime:=Set(WreathBlocksprime);
	WreathBlocksprime:={Set(v): v in WreathBlocksprime};

	return GloGenerator, h1, Block_Length,Wprime,WreathBlocksprime;

end function;


function GlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);

	EmbedImagess:=[];

	for i in [1..#alpha] do
		EmbedImagess[i]:=Evaluate(h1,alpha[i]);
	end for;

	len:=#EmbedImagess;

	Blo:={};  		//Blocks from embedding polynomial

	first:={1}; 		//First block.
	first_val:=Valuation(EmbedImagess[1]);

	Grouppedup:={1};  	//Keep track of numbers 1 thru len that have been grouped up.

	for j in [2..len] do

		if first_val eq Valuation(EmbedImagess[j]) then

			if Valuation(EmbedImagess[j] - EmbedImagess[1]) gt (first_val + 5) then

				first:=first join {j};
				Grouppedup:=Grouppedup join {j};
			end if;

		end if;

	end for;

	Blo:=Blo join {first};

	RemainIndex:={1..len} diff Grouppedup; 

	while #Grouppedup lt len do

		a:=Random(RemainIndex);
		currentblockk:={a};
		comp_val:=Valuation(EmbedImagess[a]);

		Grouppedup:=Grouppedup join {a};
		RemainIndex:=RemainIndex diff {a};

		if IsEmpty(RemainIndex) then
			alpha:=MY_LiftRoots(GloGenerator,alpha);
			return GlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
		end if;

		
		for b in RemainIndex do

			if comp_val eq Valuation(EmbedImagess[b]) then

				if Valuation(EmbedImagess[a]-EmbedImagess[b]) gt (comp_val + 5) then

					currentblockk:=currentblockk join {b};
					Grouppedup:=Grouppedup join {b};
					RemainIndex:=RemainIndex diff {b};
				end if;

			end if;

		end for;  //line 381

		if #currentblockk eq 1 then
			alpha:=MY_LiftRoots(GloGenerator,alpha);
			return GlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
		end if;

		Blo:=Blo join {currentblockk};

	end while;	

	Block_Sizes:={#a: a in Blo};

	bool:={Block_Length} eq Block_Sizes; 

	if bool then
		return Blo, alpha;

	else
		alpha:=MY_LiftRoots(GloGenerator,alpha);
                return GlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
	end if;

	return 	0;

end function;


function FindOneBlock(blo)

//Input: a set of blocks.

//Output:  the block containing 1.

        state:=false;

        while state eq false do
                b:=Random(blo);
                state:=(1 in b) or IsEmpty(blo);

                blo:=blo diff {b};
        end while;

        return b;

end function;


function ReOrderRoots(sigma,alpha)

        B:=[];

        for i in [1..#alpha] do

                B[Image(sigma,i)]:=alpha[i];

        end for;

        return B;

end function;


function MatchBlocks(WreathBlocks, SubfieldBlocks,n)

//Input:  Blocks from WreathProduct for Starting group and blocks from creating an embedding function with subfield.  The embedding function uses a global rep.  Third input is an integer n.

//Output: a permutation that send the blocks of the subfields to the blocks of the wreath product.  Permutation is a element of Sn.
	
	if #WreathBlocks ne #SubfieldBlocks then
		error "Error:  number of blocks is not equal.";
	end if;

	FirstW:=FindOneBlock(WreathBlocks);
	FirstS:=FindOneBlock(SubfieldBlocks);

	WreathBlocks:=WreathBlocks diff {FirstW};
	SubfieldBlocks:=SubfieldBlocks diff {FirstS};

	Bmap:=[];  // if 3 is sent to 4 then Bmap[3] = 4.
	
	Bmap[1]:=1;
	FirstW:=FirstW diff {1};
	FirstS:=FirstS diff {1};

	FW:=Setseq(FirstW);
	FS:=Setseq(FirstS);

	for k in [1..#FW] do
		Bmap[FS[k]]:=FW[k];
	end for;

	WreathBlocks:=Setseq(WreathBlocks);
	SubfieldBlocks:=Setseq(SubfieldBlocks);
	Sn:=Sym(n);
	
	for i in [1..#WreathBlocks] do
		W:=WreathBlocks[i];
		S:=SubfieldBlocks[i];

		//we will send the elements in S to W

		W:=Setseq(W);
		S:=Setseq(S);

		if #W ne #S then
			error "Error:  blocks aren't the same size.";
		end if;


		for j in [1..#S] do
			w:=W[j];
			s:=S[j]; 
			Bmap[s]:=w;

		end for;

	end for;

	sigma:=Sn!Bmap;

	return sigma;

end function;


function ResolventFactorization(G,H,pRootss,M,GloGenerator,p)

//Computes the relative resolvent for H<G. 
//The degrees of the irreducible factors of the resolvent are tallied and returned.

	RF:=Resolvent(G,H,pRootss,M,GloGenerator);
	FF:=pFactorDegrees(p,RF);
	FactorList:=Set(tally(FF));

	return FactorList;

end function;


function Resolvent_OrbitLengthTest(W,Candidate, H, ResolventFactorList)

	//ResolventFactorList is a set.

	G:=Candidate;
	m:=Order(W)/Order(H);

	a,b,c:=CosetAction(W,H);
	phi_G:=a(G);

	Orbb:=Orbits(phi_G);
	Orb_lengths:=[#r: r in Orbb];
	
	Orb_lengths:=Set(tally(Orb_lengths));

	Orb_lengths, TransitiveGroupIdentification(Candidate);

	return Orb_lengths eq ResolventFactorList;

end function;


function RelativeResolventFilter(W, Candidates,UpperLimit, alpha, M, GloGenerator,p)

	vprint Galois,1:"GaloisGroupMilstead: Initial Candidates", [TransitiveGroupIdentification(r): r in Candidates];

//For now, LowIndexSubgroups up to n^2

	Current:=2;
	
	while (Current lt UpperLimit) do

		vprint Galois,1:"GaloisGroupMilstead: Index", Current;

		J:=LowIndexSubgroups(W,<Current,Current>);

		for H in J do

			if NeedResolvent(W,H,Candidates) then

				FactorList:=ResolventFactorization(W,H,alpha,M,GloGenerator,p);
				Candidates:=[G: G in Candidates |Resolvent_OrbitLengthTest(W,G, H, FactorList)];
				"Candidates", [TransitiveGroupIdentification(r): r in Candidates];

				if #Candidates le 1 then

					return Candidates;
				end if;
			end if;

		end for;
		
		Current:=Current +1;		

	end while;

	return Candidates;
	
end function;



function Check_Group(G,H,alpha,M,GloGenerator,p)

//Forms the relative resolvent for the group pair H<G and checks to see if the galois group of glogenerator is contained in a conjugate of H.

	RF:=Resolvent(G,H,alpha,M,GloGenerator);
	Zp:=pAdicRing(p,900);
	Zx<x>:=PolynomialRing(Zp);

	RF:=Zx!RF;

	return HasRoot(RF);

end function;



function TestStartingGroup(W, alpha,M,GloGenerator,p)

//W is the starting group
//M is the complex bound on the pRoots alpha.

	Max_Can:=MaximalSubgroups(W); 

	C:=[r`subgroup: r in Max_Can];
	
	for c in C do
		
		bool:=Check_Group(W,c,alpha,M,GloGenerator,p);

		if bool then
			return false;
		end if;
	end for;

	return true;

end function;


// load "LastResort.magma";  //complex resolvent stuff



function CCMY_LiftRoots(GloGenerator,alpha)

	precis:=Precision(alpha[1]);
	//C:=Parent(alpha[1]);
	k:=2*precis;

	C:=ComplexField(k);
	Cx<x>:=PolynomialRing(C);
	//psi:=Cx!GloGenerator;
	psi:=GloGenerator;
	Z:=Integers();
	beta:=[];

	for j in [1..#alpha] do
		beta[j]:=HenselLift(Cx!psi,C!alpha[j],k);
		//beta[j]:=HenselLift(psi,C!alpha[j],k);
	end for;

	return beta;

end function;


function CClift_roots(Psi,alpha,k)


	//Compare precision of alpha to k.  If it's already greater than or equal to k don't lift.

	Currentprec:=Precision(Parent(alpha[1]));

	if Currentprec ge k then
		return alpha;
	end if;


	//Otherwise, lift and maintain root order

	C:=ComplexField(k);
	Cx<x>:=PolynomialRing(C);

	psi:=Cx!Psi;
	beta:=[];

	for j in [1..#alpha] do
		beta[j]:=HenselLift(Cx!psi,C!alpha[j],k);   
	end for;

        return beta;

end function;



function CCResolvent(G,H,CRoots,GloGenerator,prec)

	if (H subset G) eq false then

		error "Error: H is not a subgroup of G.";
	end if;

	if Index(G,H) eq 1 then

		error "Error:  H equals G.";
	end if;

	//Establish F, a G-relative H-invariant	

	F:=RelativeInvariant(G,H);

	// if VerifyInvariant(F,G,H) then
		vprint Galois,1:"GaloisGroupMilstead: Invariant PASSES";
	// else
		vprint Galois,1:"GaloisGroupMilstead: Invariant FAILS";
	// end if;

	CReps:=RightTransversal(G,H);
	alpha:=CRoots;

	//C:=Parent(alpha[1]);
	C:=ComplexField(prec);
	Cx<x>:=PolynomialRing(C);
	Z:=Integers();

	beta := CClift_roots(GloGenerator,alpha,prec);

	ResolventFactors:=[];

	for sigma in CReps do
		Cons:=Evaluate(F^(sigma),beta);
		Append(~ResolventFactors,x - Cons);
	end for;

	Res:=&*ResolventFactors;
	Res:=Cx!Res;
	
	coeff:=[];

	for j in [0..(Degree(Res))] do
		coeff[j+1]:=Round(Re(Coefficient(Res,j)));
	end for;

	RF:=Polynomial(Z,coeff);
	"RF",RF;

	if IsitSquareFree(RF) then
		return RF;
	end if;

	state:=false;

	while (state eq false) do
		A:=PartialTschirnhausen( GloGenerator,50);
		"Tschirn";
		L:=Parent(beta[1]);
		Ly:=PolynomialRing(L);

		A:=Ly!A;
		gamma:=[];

		for k in [1..#beta] do
			gamma[k]:=Evaluate(A,beta[k]);
		end for;

		ResolventFactors:=[];

		for sigma in CReps do
 			Cons:=Evaluate(F^(sigma),gamma);
 			Append(~ResolventFactors,x - Cons);  //Should it be x?  Shouldn't i be using Ly or is it the same as Cx?
		end for;

		Res:=&*ResolventFactors;
		Res:=Cx!Res;

		coeff:=[];

		for j in [0..(Degree(Res))] do
			coeff[j+1]:=Round(Re(Coefficient(Res,j)));
		end for;

		RF:=Polynomial(Z,coeff);
		state:=IsitSquareFree(RF);
        	
	end while;

	return RF;

end function;


function CCGlobalRepBlocks(GloGenerator,Block_Length,alpha,h1)

	EmbedImagess:=[];

	for i in [1..#alpha] do
		EmbedImagess[i]:=Evaluate(h1,alpha[i]);
	end for;

	len:=#EmbedImagess;

	Blo:={};  		//Blocks from embedding polynomial

	first:={1};

	Grouppedup:={1};  	//Keep track of numbers 1 thru len that have been grouped up.

	smallestimagee:=Min([AbsoluteValue(a): a in EmbedImagess]);

	for j in [2..len] do
		temp:=AbsoluteValue(EmbedImagess[1]-EmbedImagess[j]);
		
		if temp lt (0.001*smallestimagee) then

			first:=first join {j};
			Grouppedup:=Grouppedup join {j};
		end if;
	end for;

	Blo:=Blo join {first};

	RemainIndex:={1..len} diff Grouppedup; 

	while #Grouppedup lt len do

		a:=Random(RemainIndex);
		currentblockk:={a};

		Grouppedup:=Grouppedup join {a};
		RemainIndex:=RemainIndex diff {a};

		if IsEmpty(RemainIndex) then
			alpha:=CCMY_LiftRoots(GloGenerator,alpha);
			return CCGlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
		end if;

		for b in RemainIndex do
			temp:=AbsoluteValue(EmbedImagess[a]-EmbedImagess[b]);

			if temp lt (0.001*smallestimagee) then

				currentblockk:=currentblockk join {b};
				Grouppedup:=Grouppedup join {b};
				RemainIndex:=RemainIndex diff {b};
			end if;

		end for;

		if #currentblockk eq 1 then
			alpha:=CCMY_LiftRoots(GloGenerator,alpha);
			return CCGlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
		end if;

		Blo:=Blo join {currentblockk};

	end while;	

	Block_Sizes:={#a: a in Blo};

	bool:={Block_Length} eq Block_Sizes;

	if bool then
		return Blo, alpha;

	else
		alpha:=CCMY_LiftRoots(GloGenerator,alpha);
                return CCGlobalRepBlocks(GloGenerator,Block_Length,alpha,h1);
	end if;

	return 	0;

end function;



function CCResolventFactorization(G,H,CRoots,GloGenerator,p, prec)

//Computes the relative resolvent for H<G. 
//The degrees of the irreducible factors of the resolvent are tallied and returned.

	RF:=CCResolvent(G,H,CRoots,GloGenerator,prec);
	FF:=pFactorDegrees(p,RF);
	FactorList:=Set(tally(FF));

	return FactorList;

end function;



 
function CCRelativeResolventFilter(W, Candidates,UpperLimit, alpha, GloGenerator,p,prec)

//For now, LowIndexSubgroups up to n^2
	
	Current:=2;

	while (Current lt UpperLimit) do

		J:=LowIndexSubgroups(W,<Current,Current>);

		for H in J do

			FactorList:=CCResolventFactorization(W,H,alpha,GloGenerator,p, prec);
			Candidates:=[G: G in Candidates |Resolvent_OrbitLengthTest(W,G, H, FactorList)];
			
			if #Candidates le 1 then

				return Candidates;
			end if;

		end for;
		
		Current:=Current +1;
		
	end while;
		
	return Candidates;
	
end function;


function CCCheck_Group(G,H,alpha,GloGenerator,prec,p)

//Forms the relative resolvent for the group pair H<G and checks to see if the galois group of glogenerator is contained in a conjugate of H.

	RF:=CCResolvent(G,H,alpha,GloGenerator,prec);
	Zp:=pAdicRing(p,900);
	Zx<x>:=PolynomialRing(Zp);

	RF:=Zx!RF;

	vprint Galois,1:"GaloisGroupMilstead: RF",RF;

	return HasRoot(RF);

end function;


function CCTestStartingGroup(W, alpha,GloGenerator,prec,p)

//W is the starting group
//M is the complex bound on the pRoots alpha.

	Max_Can:=MaximalSubgroups(W); 

	C:=[r`subgroup: r in Max_Can];
	
	for c in C do
		
		bool:=CCCheck_Group(W,c,alpha,GloGenerator,prec,p);

		if bool then
			return false;
		end if;
	end for;

	return true;

end function;

// load "AMP.magma";


function TL_Auto(T,phi)

	Qp:=CoefficientRing(phi);
	n:=Degree(phi);

	if Degree(DefiningPolynomial(T)) eq 1 then
		T:=BaseRing(T);
	end if;

	if Degree(T,Qp) eq 1 then
		T:=Qp;
	end if;

	Ty<y>:=PolynomialRing(T);
	phi:=Ty!phi;	

	_, _, F:=Factorization(phi:Extensions:=true);
	TL0:=F[1]`Extension;

	if RamificationIndex(BaseRing(TL0)) gt 1 then
		V:=BaseRing(BaseRing(TL0));

		omicron:=Norm(DefiningPolynomial(TL0)); 

		if IsEisenstein(omicron) then
			omicron:=reduce_poly(omicron);
			ChangePrecision(~omicron,Precision(V));
		end if;

		A:=AutomorphismGroup(ext<V|omicron>,Qp);
	else
		A:=AutomorphismGroup(TL0,Qp);
	end if;

	return A;

end function;


function PossibleB1(Candidate,GalT1)

//Test whether there is a normal subgroup B1 of Candidate so that the quotient group is Isomorphic to Gal(T1/K).  If so, returns possible B1.

	H:=Candidate;
	
	ord:=Order(GalT1);
	Pos:=LowIndexSubgroups(H,<ord,ord>);

	if IsEmpty(Pos) then
		return false, [0];
	end if;
	
	B:=[];  //holds the possible subgroups B1 if any.

	for BB in Pos do
		bool:=IsIsomorphic(quo<H|BB>,GalT1) and IsNormal(H,BB);
		if bool then
			Append(~B,BB);
		end if;
	end for;

	if #B gt 0 then
		return true, B;
	else
		return false, [0];
	end if;

	return 0;

end function;



function PossibleD(Candidate,DL1,Gal1)

//Test whether: 1. there is a subgroup D1 of Candidate so that [Candidate:D1] = [L1:K].
//		2. there is a normal subgroup D0 of D1 so that D1/D0 is iso to Gal(L0,L1) = Gal1


//If true, output set {<D1,[D0]>} where [D0] is a list.

//Notation:  DL1 is Degree(L1,K).


	H:=Candidate;

	Ord:=DL1;

	D1:=LowIndexSubgroups(H,<Ord,Ord>);  //Possible D1

	if IsEmpty(D1) then
		return false, [0];
	end if;
	
	D_tuples:={};

	for d1 in D1 do

		Pos:=LowIndexSubgroups(d1,<Order(Gal1),Order(Gal1)>);
		d0list:=[];
		for d0 in Pos do

			bool:=IsIsomorphic(Gal1, quo<d1|d0>);

			if bool then
				Append(~d0list,d0);
			end if;

		end for;

		if #d0list gt 0 then 
			D_tuples:=D_tuples join {<d1,d0list>};
		end if;

	end for;


	if IsEmpty(D_tuples) then
		return false, [0];
	end if;

	return true, D_tuples;


end function; 




function AMP_Max_Tame_Test(Candidate, TameGrp, p, B1)

	PP:=TameGrp;
	
	if Order(TameGrp) eq 1 then
		state:= IspGroup(Candidate,p);

		if state eq false then
			return false, [0];
		end if;

	end if;	


	ord:=Order(PP);

	Pos:=LowIndexSubgroups(Candidate,<ord,ord>);

	if IsEmpty(Pos) then 
		return false, [0];
	end if;

	B:=[];  //holds the possible subgroups S if any.

	for S in Pos do 
		bool:=IsIsomorphic(PP, quo<Candidate|S>) and IsNormal(Candidate,S) and IspGroup(S,p) and (S subset B1);

		if bool then

			if IsNormal(B1,S) then
				Append(~B,S);
			end if;
		end if;
	end for;

	if #B gt 0 then
		return true, B;
	else
		return false, [];
	end if;

	return 0;

end function;



function PrevNormalTest(Candidate,N1,B1)

//Let N1 be Gal(L1/K), Candidate a possible Gal Grp.
//Test whether there is a normal subgroup U of Candidate so that Candidate/U is isomorphic to N1.  Also make sure U is a normal subgroup of B1

//Output:  first a boolean value for whether the candidate passes.  If it passes (true) then the second output is the list of possible U.

	H:=Candidate;
	ord:=Order(N1);
	Pos:=LowIndexSubgroups(H,<ord,ord>);

	if IsEmpty(Pos) then
		return false, [0];
	end if;

	B:=[];  //holds the possible subgroups U if any.

	for U in Pos do
		bool:=IsIsomorphic(quo<H|U>,N1 ) and IsNormal(H,U) and (U subset B1);

		if bool then

			if IsNormal(B1,U) then
				Append(~B,U);
			end if;
		end if;
	end for;

	if #B gt 0 then
		return true, B;
	else
		return false, [0];
	end if;

	return 0;

end function;


function PossibleCandidateOrders(N1,s1,T,T1,DL1)
	//DL1 is Degree(L1,K).

        Ord:=[];   //Possible Group orders.  The output of the function

//First find possible upper bound on w.
        K:=PrimeRing(T1);
vprint Galois,1:"GaloisGroupMilstead: [N1:K]",Order(N1);
vprint Galois,1:"GaloisGroupMilstead: [T1:K]",Degree(T1,K);
vprint Galois,1:"GaloisGroupMilstead: [T:K]",Degree(T,K);
vprint Galois,1:"GaloisGroupMilstead: [T:T1]",Degree(T,K)/Degree(T1,K);
	p:=Prime(K);
	temp:=Lcm(Order(N1),Degree(T,K));
        tamemult:=Integers()!(Order(N1)/Degree(T1,K))*Integers()!(Degree(T,K)/Degree(T1,K));
        wildbase:=p;
        stemtame := DL1/p^Valuation(DL1,p);
        wildexps:=[0..Integers()!(Order(N1)/Degree(T1,K))*s1*stemtame];
vprint Galois,1:"GaloisGroupMilstead: wildexps",wildexps;
	for w in wildexps do
          for t in Divisors(tamemult) do
		finall:=temp * wildbase^w *t;
		finall:=Integers()!finall;
		Append(~Ord,finall);
          end for;
	end for;
vprint Galois,1:"GaloisGroupMilstead: Possible Orders",Ord;
	return Ord;

end function;


function TwoSegmentStartingGroup(L1,L0,Prev_Grp)

//Prev_Grp is Gal(L1/K) if this isn't the first iteration.  Otherwise it is the trivial group.

//Output:  W=Starting Group as a wreath product
	 //N1 = Gal(L1/K)
	 //WreathBlocks
	//Gal1 is Gal(L0/L1)

	if Order(Prev_Grp) gt 1 then
		N1:=Prev_Grp;
	else
		if IsWildlyRamified(L1) then
			N1:=One_Segment(DefiningPolynomial(L1));
		else
			N1:=TameGaloisGroup(L1);
		end if;
	end if;

	if IsWildlyRamified(L0) then
		Gal1:=One_Segment(DefiningPolynomial(L0));
	else
		Gal1:=TameGaloisGroup(L0);
	end if;


	W:=WreathProduct(Gal1,N1);
	
	ba, im, ker:=BlocksAction(W,{1..Degree(Gal1)});
	WreathBlocks:=Orbits(ker);

	WreathBlocks:=Set(WreathBlocks);
	WreathBlocks:={Set(w): w in WreathBlocks};

	return W, N1, WreathBlocks,Gal1;


end function;


function Two_Segments_Filtering(phi,L1, L0, Prev_Grp, T,T1)
vprint Galois,1:"GaloisGroupMilstead: TwoSegmentsFiltering";

	Num_Remaining:=[];
	
	Qp:=CoefficientRing(phi);
	p:=Prime(Qp);

	n:=Degree(phi);
	K:=CoefficientRing(phi);

	W, N1, WreathBlocks,Gal1:=TwoSegmentStartingGroup(L1,L0,Prev_Grp); 

	//TransitiveGroupIdentification(Gal1), "Gal1";

	s1:=Valuation(Degree(DefiningPolynomial(L0)),p);
	DL1:=Degree(L1,K);
	Pos_Order:=PossibleCandidateOrders(N1,s1,T,T1,DL1);

	if n lt 35 then
vprint Galois,1:"GaloisGroupMilstead: Finding transitive subgroups of", W;
		Can:=Subgroups(W:IsTransitive:=true);
		C:=[r`subgroup: r in Can];
		InitialCandidates:=CandidateClasses(C);  //associative array

		First:=Keys(InitialCandidates);
		Num_Remaining[1]:=#First;

		Second:=[c: c in First| Order(c) in Pos_Order];
		Num_Remaining[2]:=#Second;
	else 
		min_order:=Min(Pos_Order);
		max_order:=Max(Pos_Order); //Use list of possible orders to set parameters for calculating initial group candidates as subgroups of W.
vprint Galois,1:"GaloisGroupMilstead: Finding transitive subgroups of", W,"of order dividing",max_order;

		max_order:=Min(max_order,Order(W));  //Need to limit ourselves to size of W.

		Can:=Subgroups(W:IsTransitive:=true, OrderMultipleOf:=min_order, OrderDividing:=max_order );
		C:=[r`subgroup: r in Can];
		InitialCandidates:=CandidateClasses(C);  //associative array

		First:=Keys(InitialCandidates);
		Num_Remaining[1]:=#First;

		Second:=First;
		Num_Remaining[2]:=#Second;
	end if;
vprint Galois,1:"GaloisGroupMilstead: #Remaining 2:",Num_Remaining;

	PParity:=PolynomialParity(phi);
	third:=[G: G in Second | (PParity eq GroupParity(G))];
	
	Num_Remaining[3]:=#third;

vprint Galois,1:"GaloisGroupMilstead: #Remaining 3:",Num_Remaining;
	if #third eq 1 then
		return Num_Remaining, InitialCandidates, third, W, WreathBlocks;
	end if;
	
	Autoo:=AutomorphismGroup(ext<Qp|phi>,Qp);
	Fourth:=[];

 	if Order(Autoo) eq n then
		Num_Remaining[4]:=1;
		return Num_Remaining, InitialCandidates, [Autoo], W, WreathBlocks;
	end if;

	for i in [1..#third] do
		G:=third[i];
		Cen:=Centralizer(Sym(n),G);
	
		if IsIsomorphic(Autoo,Cen) then
			Append(~Fourth,G);
		end if;

	end for;

	Num_Remaining[4]:=#Fourth;
vprint Galois,1:"GaloisGroupMilstead: #Remaining 4:",Num_Remaining;

	if #Fourth eq 1 then
		return Num_Remaining, InitialCandidates, Fourth, W, WreathBlocks;
	end if;


	BTuples:=[];

	GalT1:=TameGaloisGroup(T1);

	for n in Fourth do
		bool,B1:=PossibleB1(n,GalT1);
		
		if bool then
			Append(~BTuples,<n,B1>);
		end if;
	end for;

	TameGrp:=TameGaloisGroup(T);
	Fifth:=[];  //after using in Num_Remaining, rename it BSTuples

	for m in BTuples do  //now over BTuples
		
		Candidate:=m[1];

		B:=m[2];  		//Possible B's.  Loop over this.

		BS_set:={};  //Set of tuples <B1,[S]>
		for B1 in B do

			bool,S:=AMP_Max_Tame_Test(Candidate, TameGrp, p, B1);
			//bool;
	
			if bool then
				BS_set:=BS_set join {<B1,S>};
			end if;
		end for;
		
		if #BS_set gt 0 then
	
			Append(~Fifth,<Candidate,BS_set>);
		end if;

	end for;

	Num_Remaining[5]:=#Fifth;

vprint Galois,1:"GaloisGroupMilstead: #Remaining 5:",Num_Remaining;
	if #Fifth eq 1 then
		return Num_Remaining, InitialCandidates, [t[1]: t in Fifth], W, WreathBlocks;
	end if;
	
	BSTuples:=Fifth;

	//return Num_Remaining;

	Sixth:=[];  	//after using in Num_Remaining, rename it BSU_Tuples.  Contains B1, S, and U.

	for m in BSTuples do
		Candidate:=m[1];
		CurrentSet:=m[2];  //of the form {<B1, [S]>}.

		BSU_set:={};

		for C in CurrentSet do
			B1:=C[1];
			Possible_S:=C[2];

			bool,U:=PrevNormalTest(Candidate,N1,B1);

			if bool then
				BSU_set:= BSU_set join {<B1,Possible_S,U> };
			end if;	

		end for;

		if #BSU_set gt 0 then
			Append(~Sixth,<Candidate,BSU_set>);
		end if;

	end for;

	Num_Remaining[6]:=#Sixth;

vprint Galois,1:"GaloisGroupMilstead: #Remaining 6:",Num_Remaining;
	if #Sixth eq 1 then
		return Num_Remaining, InitialCandidates, [s[1]: s in Sixth], W, WreathBlocks;
	end if;
	
	BSU_Tuples:=Sixth;

	//Change the next four lines

	Seventh:=[];

	for t in BSU_Tuples do

		H:=t[1];
		CurrentSet:=t[2];

		New_set:={};   // {<B1, SU_set>}

		for C in CurrentSet do
			B1:=C[1];
			Possible_S:=C[2];
			Possible_U:=C[3];
			SU_set:={};			//{<s,u>}

			for s in Possible_S do

				for u in Possible_U do
					bool:=IsIsomorphic(quo<B1|s>, quo<u|(s meet u)>);
				

					if bool then
						SU_set:=SU_set join {<s,u>};
					end if;

				end for;

			end for;

			if #SU_set gt 0 then
				New_set:=New_set join {<B1, SU_set>};

			end if;


		end for;

		if #New_set gt 0 then
			Append(~Seventh,<H,New_set>);
		end if;

	end for;

	Num_Remaining[7]:=#Seventh;

vprint Galois,1:"GaloisGroupMilstead: #Remaining 7:",Num_Remaining;
	if #Seventh eq 1 then
		return Num_Remaining, InitialCandidates, [s[1]: s in Seventh], W, WreathBlocks;
	end if;
	
	Eighth:=[];

	for v in Seventh do

		H:=v[1];
		Btupleset:=v[2];

		bbool, D:=PossibleD(H,DL1,Gal1);

		if bbool then
			Append(~Eighth,<H,Btupleset,D>);
		end if;

	end for;	

	Num_Remaining[8]:=#Eighth;
	
vprint Galois,1:"GaloisGroupMilstead: #Remaining 8:",Num_Remaining;
	if #Eighth eq 1 then
		return Num_Remaining, InitialCandidates, [s[1]: s in Eighth], W, WreathBlocks;
	end if;
	
	Ninth:=[];

	for r in Eighth do

		H:=r[1];

		BB:={};   //will be {<B1,s,u,D1,D0>}

		for T in r[2] do
			B1:=T[1];
			su_pairs:=T[2];

			for w in su_pairs do
				s:=w[1];
				u:=w[2];

				for z in r[3] do
		        		D1:=z[1];
                        		if u subset D1 then
						for D0 in z[2] do
							Group1:=quo<(s meet D1)|(s meet D0)>;
							bool1:= (Order(Group1) eq p^s1) and IsElementaryAbelian(Group1);
							Group2:=quo<(s meet u)|((s meet u) meet D0)>;
							bool2:=(Order(Group2) le p^s1) and IsElementaryAbelian(Group2) and IspGroup(Group2,p);

							TN1L0_TN1 := #(u meet s)/#(D0 meet u meet s);
							TN1_K := #H/#(u meet s);
							e0 := DL1/p^Valuation(DL1,p);
							N1_T1 := Order(N1)/Degree(T1,K);
							possible_orders := {TN1_K*p^w:w in [Valuation(TN1L0_TN1,p)..e0*N1_T1*Valuation(TN1L0_TN1,p)]};

                                        		bool3:=#H in possible_orders;
										
							if (bool1 and bool2 and bool3) then
								BB:=BB join {<B1,s,u,D1,D0>};
							end if;

			 			end for;

                            		end if;
				end for;

			end for;
		end for;

		if #BB gt 0 then
			Append(~Ninth,<H,BB>);
		end if;
	end for;

	Num_Remaining[9]:=#Ninth;
	 
	Candidates:=[Ninth[i][1] : i in [1..#Ninth]];     //Later: after 8th criteria.

vprint Galois,1:"GaloisGroupMilstead: #Remaining 9:",Num_Remaining;
	if #Candidates eq 1 then
		return Num_Remaining, InitialCandidates, Candidates, W, WreathBlocks;
	end if;

	Candidates:=GlobalGalFilter(Candidates,phi);

	Num_Remaining[10]:=#Candidates;

vprint Galois,1:"GaloisGroupMilstead: #Remaining 10:",Num_Remaining;
	if #Candidates eq 1 then
		return Num_Remaining, InitialCandidates, Candidates, W, WreathBlocks;
	end if;

	Sn:=Sym(n);

	if NeedResolvent(Sn,DirectProduct(Sym(2),Sym(n-2)),Candidates) then
		tau:=phi;

		psi:=S2Resolvent(tau);

		state:=IsitSquareFree(psi);

		while (state eq false) do
			tau:=Tschirnhausen(tau);
			psi:=S2Resolvent(tau);
			state:=IsitSquareFree(psi);

		end while;

		FF:=pFactorDegrees(p,psi);
		S2ResolventFactorList:=Set(tally(FF));

		New_Candidates:=[G: G in Candidates |S2_OrbitLengthTest(G,tau,S2ResolventFactorList)];
				
	else

		New_Candidates:=Candidates;   //keep same list since not affecting it.
	end if;

	Num_Remaining[11]:=#New_Candidates;		
	
vprint Galois,1:"GaloisGroupMilstead: #Remaining 11:",Num_Remaining;
	if #New_Candidates eq 1 then
		return Num_Remaining, InitialCandidates, New_Candidates, W, WreathBlocks;
	end if;



	//time Autoo:=TL_Auto(T,phi);

	//Eleventh:=[];

	//for y in Ninth do

		//H:=y[1];

		//if H in New_Candidates then
			//Poss:=[];

			//tuppSet:=y[2];

			//for v in tuppSet do
				//s:=v[2];
				//D0:=v[5];

				//LBool:=IsIsomorphic(Autoo, quo<Normalizer(H,(s meet D0))|(s meet D0)>);
				//Append(~Poss,LBool);

			//end for;

			//if true in Poss then
				//Append(~Eleventh,H);
			//end if;

		//end if;

	//end for;

	//Eleventh:=Set(Eleventh);
	//Num_Remaining[11]:=#Eleventh;
	
	//return Num_Remaining, InitialCandidates, Eleventh, W, WreathBlocks;

vprint Galois,1:"GaloisGroupMilstead: Numbers:",Num_Remaining; //SP
vprint Galois,1:"GaloisGroupMilstead: Orders",[Order(g):g in Candidates];

	if NeedResolvent(Sn,DirectProduct(Sym(1),DirectProduct(Sym(1),Sym(n-2))),New_Candidates) then  		//H:=S_1 x S_1 x S_{n-2}

		lrResolventFactorList:=lrFact(phi);
		C1:=[G: G in New_Candidates| lr_OrbitLengthTest(G,phi,lrResolventFactorList)];
		
	else
		C1:=New_Candidates;
		
	end if;

	Num_Remaining[#Num_Remaining +1]:=#C1;

vprint Galois,1:"GaloisGroupMilstead: #Remaining C1:",Num_Remaining;
//a:=0; 1/a;
	if #C1 eq 1 then
		return Num_Remaining, InitialCandidates, C1, W, WreathBlocks;
	end if;


	if NeedResolvent(Sn,DirectProduct(Sym(3),Sym(n-3)),C1) then 				//H:=S_3 x S_{n-3}

		S3ResolventFactorList:=tpFact(phi);
		C2:=[G: G in C1 |S3_OrbitLengthTest(G,phi,S3ResolventFactorList) ];

	else
		C2:=C1;
	end if;

	Num_Remaining[#Num_Remaining +1]:=#C2;

vprint Galois,1:"GaloisGroupMilstead: #Remaining C2:",Num_Remaining;
	if #C2 eq 1 then
		return Num_Remaining, InitialCandidates,C2, W, WreathBlocks;
	end if;



	if NeedResolvent(Sn, DirectProduct(Sym(2),DirectProduct(Sym(1),Sym(n-3))), C2) then
		LRResolventFactorList:=LRFact(phi);
		C3:=[G: G in C2 | LR_OrbitLengthTest(G,phi,LRResolventFactorList) ];

	else

		C3:=C2;
	end if;

	Num_Remaining[#Num_Remaining +1]:=#C3;


vprint Galois,1:"GaloisGroupMilstead: #Remaining C3:",Num_Remaining;
	return Num_Remaining, InitialCandidates, C3, W, WreathBlocks;

end function;



function Two_Segments_Eisenstein(phi,L1, L0, Prev_Grp, T,T1)
vprint Galois,1:"GaloisGroupMilstead: Two Segments";
//First, we reduce coefficients of phi.
	n:=Degree(phi);
	p:=Prime(PrimeRing(T));

	phi:=reduce_poly(phi);	
	Qp:=CoefficientRing(phi);	
	ChangePrecision(~phi,Precision(Qp));

//Second, determine the Candidates and filter.

	Num_Remaining, InitialCandidates, Candidates, W, WreathBlocks:=Two_Segments_Filtering(phi,L1, L0, Prev_Grp, T,T1); //InitialCandidates is an array.  Candidates is the collection of remaining keys.

vprint Galois,1:"GaloisGroupMilstead: Numbers:",Num_Remaining; //SP
vprint Galois,1:"GaloisGroupMilstead: Orders",[Order(g):g in Candidates];
// SP	//if (n le TransitiveGroupDatabaseLimit()) and (#Candidates eq 1) then
	if (#Candidates eq 1) then
		//Num_Remaining;
		//phi;
		//TransitiveGroupIdentification(Candidates[1]);
vprint Galois,1:"GaloisGroupMilstead: #Remaining",#Num_Remaining,Num_Remaining;
		return Candidates[1],Num_Remaining;						//Early exit
	end if;
	if n gt TransitiveGroupDatabaseLimit() then 
		Candidates:=Flat([Setseq(InitialCandidates[a]): a in Candidates]);  //Put all groups in remaining isomorphism classes in a list.
	end if;
// SP  
//print Candidates;
//print "bis hier und nicht weiter";
//a:=0; 2/a;
vprint Galois,1:"GaloisGroupMilstead:  ----------- weiter ";
//More than one Candidate remains so we get the necessary data to form resolvents.

	GloGenerator, h1, Block_Length,Wprime,WreathBlocksprime:=GGlobalRepresentation(phi,L0);   //h1 is embedding map.
	alpha:=pAdicRoot_Approx(GloGenerator,1000);

	SubfieldBlocks, alpha:=GlobalRepBlocks(GloGenerator, Block_Length,alpha,h1);

	M:=FujiwaraBound(GloGenerator);

	//Here we need to compare W and Wprime since the rest depends on the correct block structure.

	if MyIsConjugate(W,Wprime) then                       //Check if "same" group. Maybe change so you check equal.

		sigma:=MatchBlocks(WreathBlocks,SubfieldBlocks,n);
		alpha:=ReOrderRoots(sigma,alpha);	

		if W in Candidates then			//If necessary check the starting group.

			if  TestStartingGroup(W,alpha,M,GloGenerator,p) then
				return W,Num_Remaining;
			else 
				Exclude(~Candidates,W);
			end if;
		end if;

		if (n le TransitiveGroupDatabaseLimit()) and (#Candidates eq 1) then
			return Candidates[1],Num_Remaining;
		end if;

		//Try and narrow down the candidates by looking at resolvents for low index subgroups of W.

		UpperLimit:=Max([Index(W,c):c in Candidates]) +1;  //just the max apparently isn't enough
		//UpperLimit:=n^2;
	
		Candidates:=RelativeResolventFilter(W, Candidates,UpperLimit, alpha, M, GloGenerator,p);

		if (n le TransitiveGroupDatabaseLimit()) and (#Candidates eq 1) then
			return Candidates[1],Num_Remaining;
		end if;

		return Candidates, Num_Remaining;

	end if;

	//Case where W is a subset of Wprime

	if (W subset Wprime) then

		vprint Galois,1:"GaloisGroupMilstead: BURP!!! ";

		sigma:=MatchBlocks(WreathBlocksprime,SubfieldBlocks,n);
		alpha:=ReOrderRoots(sigma,alpha);

		UpperLimit:=Max([Index(Wprime,c):c in Candidates]) +1;  //just the max apparently isn't enough
		//UpperLimit:=n^2;
	
		Candidates:=RelativeResolventFilter(Wprime, Candidates,UpperLimit, alpha, M, GloGenerator,p);

		if (n le TransitiveGroupDatabaseLimit()) and (#Candidates eq 1) then
			return Candidates[1],Num_Remaining;
		end if;

		return Candidates, Num_Remaining;

	end if;


	//Case where W is not a subset of Wprime directly.

	error "This Case has not been implemented yet.";

	return Candidates, Num_Remaining;

end function;


function Three_Segments_Eisenstein(L)
vprint Galois,1:"GaloisGroupMilstead: Three Segments";

	Qp:=PrimeRing(L[1]);

	if #L ne 3 then
		error "Error: The number of segments is wrong";
	end if;

	if IsTamelyRamified(L[1]) then   //Same as seeing if there is a horizontal segment.

			T1:=L[1];   //same in Pauli's code.
			T:=Max_Tame_Subextension(DefiningPolynomial(L[2],Qp));
	else
			T1:=Max_Tame_Subextension(DefiningPolynomial(L[1],Qp));
			T:=Max_Tame_Subextension(DefiningPolynomial(L[2],Qp));
	end if;

	N1:= Two_Segments_Eisenstein(DefiningPolynomial(L[2],Qp),L[1], L[2], CyclicGroup(1), T,T1);	//will be prevgrp after first

	//Now for third segment

	TT:=Max_Tame_Subextension(DefiningPolynomial(L[3],Qp));

	return Two_Segments_Eisenstein(DefiningPolynomial(L[3],Qp),L[2], L[3], N1, TT,T);

end function;




function Eisenstein_Case(phi:global:=false)
vprint Galois,1:"GaloisGroupMilstead: Eisenstein_Case";
	if IsEisenstein(phi) eq false then
		error "Error:  Input polynomial is not Eisenstein";
	end if;
	
	phi:=reduce_poly(phi);
	n:=Degree(phi);

	Qp:=CoefficientRing(phi);
	p:=Prime(Qp);

	prec:=Precision(Qp);

	ChangePrecision(~phi,prec);

	//Catch precision if it's too low for discriminant.

	ramification_polygon,rho:=My_RamificationPolygon(phi);
	vertices, slopes:=Recover_Vertices_Slopes(ramification_polygon);

	my_discriminant:=RecoverDiscriminant(vertices);		//Find valuation of discriminant.
	my_discriminant:=Integers()!my_discriminant;

	if prec le my_discriminant then
vprint Galois,1:"GaloisGroupMilstead: Precision incresed to",p,"^",my_discriminant + 3;
		
		Zp:=pAdicRing(p, my_discriminant + 3);
		Zx:=PolynomialRing(Zp);
		return Eisenstein_Case(Zx!phi);
	end if;
	
	//Early Exits.

	if Valuation(n,p) eq 0 then

		return TameGaloisGroup(StemField(phi)); 
	end if;


	if n eq p then
		phi:=MakeEisenstein(phi);
		return One_Segment(phi);
		
	end if;

	K:=ext<Qp|phi>;		//stem field

	if (Number_Of_Segments(phi) eq 1) and (n eq (p^(Valuation(n,p)))) then
		 return One_Segment(phi);
	end if;
vprint Galois,1:"GaloisGroupMilstead: AutomorphismGroup";
	Autoo:=AutomorphismGroup(K,Qp);  //Maybe not do this here.

	if Order(Autoo) eq n then
		return Autoo;  //We are in the normal case here.
		
	end if;

	L:=Ramification_Segments(K);

	if #L eq 2 then
		L1:=L[1];
		L0:= L[2];

		if IsTamelyRamified(L[1]) then   //Same as seeing if there is a horizontal segment.

			T1:=L1;   //same in Pauli's code.
			T:=Max_Tame_Subextension(DefiningPolynomial(L0,Qp));
		else
			T1:=Max_Tame_Subextension(DefiningPolynomial(L1,Qp));
			T:=Max_Tame_Subextension(DefiningPolynomial(L0,Qp));
		end if;

		Prev_Grp:=CyclicGroup(1);

                G, L := Two_Segments_Eisenstein(DefiningPolynomial(L0,Qp),L1, L0, Prev_Grp, T,T1);
		return G;
	end if;

	if #L eq 3 then

		return Three_Segments_Eisenstein(L);

	end if;
error "GaloisGroupMilstead: not implemeneted yet";
	return 0;

end function;

// load "Sinclair.magma";


function ProcessVertices(lis)

//Input: List of vertices.

//Purpose: remove the vertices that aren't endpoints.

	N:=NewtonPolygon(lis);	
	vertices:=LowerVertices(N);

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];	//get rid of point with x-coordinate=0. 
	end if;

	return vertices;

end function;



function MySlope(tup1,tup2)

	//Input:  tup1=(x1,y1), tup2=(x2,y2).

	//Output: slope of line through the two points.

	x1:=tup1[1];
	y1:=tup1[2];

	x2:=tup2[1];
	y2:=tup2[2];

	if x1 eq x2 then
		return "Undefined slope";  //Type is MonStgElt
	end if;

	m:=(y2-y1)/(x2-x1);

	return Rationals()!m;

end function;



function FindVertices(slopelist, segmentlengths)

	//Given the slopes and segment lengths find the vertices of the ramification polygon.

	n:=1 + &+segmentlengths;

	vertices:=[<n,0>];  //last vertex.  We will work backwards and then reverse this list.

	segmentlengths:=Reverse(segmentlengths);
	slopes:=Reverse(slopelist);

	Currentx:=n;
	Currenty:=0;   //Now while loop
	i:=1;

	while i lt (#segmentlengths +1) do

		Newx:=Currentx-Integers()!(segmentlengths[i]);
		m:=slopes[i];		

		Newy:= Currenty + Integers()!(m*(Newx - Currentx));

		Append(~vertices,<Newx, Newy>);

		Currentx:=Newx;
		Currenty:=Newy;

		i:=i+1;

	end while;
	
	vertices:=Reverse(vertices);

	return vertices;

end function;



function NextLowerPolygon(vertexlist, slopelist)


//Let L > L1 > ...>Ll >K be the tower of subfields from the ramification polygon
//of L/K.  This function determines the ramification polygon, in vertices and slopes, of L1/K.

//Input:  the vertices and slopes of the ramification polygon of L/K.

//Output: the vertices and slopes of the ramification polygon of L1/K.  Also the shrinkfactor ps1.


	if #slopelist eq 1 then
		error "Error: Polygon only has 1 segment.  There is no subfield chain.";
	end if;


//First, make sure the inputted vertex list is in the right form.

	if vertexlist[1][1] eq 0 then
		vertexlist:=vertexlist[2..#vertexlist];		//get rid of point with x-coordinate=0
		slopelist:=slopelist[2..#slopelist];		//get rid of first segment.
	end if;

	if vertexlist[1][1] ne 1 then
		error "Error: polygon doesn't have a vertex with x-coordinate 1.";
	end if;  //Fine up to here.

	vertexlist:=vertexlist[2..#vertexlist];  
	slopelist:=slopelist[2..#slopelist];	//get rid of first segment.

	ps1:=vertexlist[1][1];  	//first endpoint in new polygon


//Now determine the lengths of the new segments.

	NewLengths:=[];

	oldxcoords:=[Integers()!vertexlist[i][1]: i in [1..#vertexlist]];

	for i in [1..#slopelist] do
		c:=oldxcoords[i+1] - oldxcoords[i];
				
		NewLengths[i]:=Integers()!(c/ps1);
	end for;

	vertices:=FindVertices(slopelist, NewLengths);
		
	return  vertices, slopelist,ps1;


end function;



function  Possible_Phi0(K,n)

	//K is base field, n is degree of polynomial.  n can be recovered from polygon.

	p:=Prime(PrimeRing(K));

	m:=Valuation(n,p);

	if (n eq p^m) then
		return {p};
	else
		_K_,vk:=ResidueClassField(K);	//vk is Mapping from: RngPad: K to FldFin: _K_

		Kx,vkx:=UnitGroup(_K_);		//vkx is Mapping from: GrpAb: Kx to FldFin: _K_

		if IsCyclic(Kx) then

			a:=Setseq(Generators(Kx))[1];
			den:=sub<Kx|n*a>;

			S,vS:=quo<Kx |den>;  	//vS is Mapping from: GrpAb: Kx to GrpAb: S.  

			//Correct to here.

			Delta:={(vkx(s@@vS))@@vk  :s in S};   //Map each element of S to K (through _K_ and Kx)

			return {delta*p: delta in Delta};
		else
			
			error "Error: Something has gone terribly wrong.";
		end if;
		
	end if;		

end function;


function Ram_OneSegment(Segmentslope, residual_polynomial,BaseField,PolyDegree, leftend)

	if Type(residual_polynomial) eq SeqEnum then
		error "Error: Second input should be a polynomial.";
	end if;

	K:=BaseField;
	n:=PolyDegree;
	A:=residual_polynomial;

	if Type(K) eq FldPad then

		error "Coefficient Ring must be a ring, not a field.";
	end if;
	
	Fp:=ResidueClassField(PrimeRing(K));
	p:=#Fp;
	m:=Valuation(n,p);

	Fq:=ResidueClassField(K);
	q:=#Fq;
	Fx:=PolynomialRing(Fq);

	e:=Denominator(-Segmentslope);
	h:=Numerator(-Segmentslope);

	A:=Fx!A;
	f1:=Segmental_Inertia(A);
	f:=LCM(f1,UnityIndex(K,e));

	b, b_tilde:=Difference_XGCD(h,e);
	a, a_tilde:=Difference_XGCD(e,n);

	Fqf:=ResidueClassField(UnramifiedExtension(K,f));
	Fqf_unitgroup, unitmap:= UnitGroup(Fqf);
	zeta:=unitmap(Fqf_unitgroup.1);

	Fz<z>:=PolynomialRing(Fqf);
	A:=Fz!A;
	
	//Now to find the roots of A and r.

	RR:=Roots(A);
	u:=[r[1]^leftend: r in RR];
	u1:=u[1];

	if Log(u1) eq 0 then
		r_prime:=0;
	else
		r_prime:=Integers()!(Log(u1^b)/Log(zeta));
	end if;

	
	if r_prime in  Seqset([0..(e-1)]) then
		r:=r_prime;

	else
		r:=(r_prime mod e);

	end if;	


	//Finding M

	Fqf_plus,ADDMap:=AdditiveGroup(Fqf); //ADDMap maps Fqf_plus to Fqf
	M:=sub<Fqf_plus|Identity(Fqf_plus)>;
	i:=1;

	state:=Order(M) eq p^m;


	while state eq false do

		Cons:=u[i]/(zeta^(r*h));
		Cons:=Fqf!Cons;
		psi:=z^e -Cons;

		RRi:=Roots(psi);
		ui:=[r[1]: r in RRi];
		gens:=[(a*s)@@ADDMap: s in ui];
		M:=sub<Fqf_plus|gens cat Setseq(Generators(M))>;
		i:=i+1;

		state:=Order(M) eq p^m;

	end while;
	
	MGens:=Generators(M);
	V:=VectorSpace(Fp,#MGens);
	B:=Basis(sub<V|[V!Eltseq(b): b in MGens]>);

	k:=(r*(q-1))/e;
	l:=(-1+q^f)/e;

	k:=Integers()!k;
	l:=Integers()!l;

	Expss:=[];	

	for j in [1..Dimension(V)] do

		current:=M.j;
		current:=ADDMap(current);

		if Log(current) eq 0 then
			Expss[j]:=0;
		else
			Expss[j]:=Integers()!(Log(current)/Log(zeta));

		end if;				

	end for;


	S_rows:=[];
	T_rows:=[];


	for s in [1..Dimension(V)] do		//Finding the rows of S using the automorphism s defined on powers of zeta

		j:=Expss[s];

		temp:=zeta^(l*h +j);
		temp:=temp@@ADDMap;
		temp:=M!temp;
		temp:=V!Eltseq(temp);

		temp:=ElementToSequence(temp);

		S_rows[s]:=temp;

	end for;


	for t in [1..Dimension(V)] do		//Find the rows of T using the automorphism t defined on powers of zeta

		j:=Expss[t];

		temp:=zeta^(h*k +q*j);
		temp:=temp@@ADDMap;
		temp:=M!temp;
		temp:=V!Eltseq(temp);

		temp:=ElementToSequence(temp);

		T_rows[t]:=temp;  		

	end for;


	S:=Matrix(GF(p), S_rows);
	T:=Matrix(GF(p), T_rows);

	GGLL:=GL(m,Fp);

	S:=GGLL!S;
	T:=GGLL!T;

	Matgenn:=sub<GGLL|S,T>;

	G:=AffineGroup(Matgenn);
	
	return G;

end function;




function Possible_Auto_Orders(vertexlist, slopelist, residual_polynomials)

	if Type(residual_polynomials) eq Tup then

		residual_polynomials:=TupleToSequence(residual_polynomials);
	end if;

	if Type(residual_polynomials) ne SeqEnum then
		error "Error: Third input should be a list.";
	end if;

	Count:={1};  //1 root from infinite segment, which is typically tossed.

	vertices:=vertexlist;
	slopes:=slopelist;
	
	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;


	a:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices
	
	for i in [1..#slopes] do

		if (slopes[i] in Integers()) then

			di:=a[i+1]-a[i];    //length of segment
			Ai:=residual_polynomials[i];

			if HasRoot(Ai) then
				
				RRoots:=Roots(Ai);

				for j in [1..#RRoots] do

					r:=RRoots[j][2];  //multiplicity of j-th root.

					if r eq 1 then
						Count:={a +1: a in Count};  //if root is simple add 1 to list
					else
						Count:={a + b: a in Count, b in [0..r] }; // otherwise add between 0 and r to list.

					end if;

				end for;
			end if;


		end if;

	end for;

	n:= Integers()!Max(a); //Degree of polynomial

	//ALL THAT'S LEFT IS TO THROW OUT ORDERS THAT DON'T DIVIDE n=DEG(POLY)

	return Set([b: b in Count|IsDivisibleBy(n,b)]);

end function;



function Possible_AutoGroups(vertexlist, residual_polynomials)

	vertices:=ProcessVertices(vertexlist);

	slopes:=[];

	for i in [1..(#vertices -1)] do  	//Find slopes using the function MySlope

		tup1:=vertices[i];
		tup2:=vertices[i+1];
		s:=MySlope(tup1,tup2);

		Append(~slopes,s);
	end for;


	ords:=Possible_Auto_Orders(vertices, slopes, residual_polynomials);

	possi:=[];	

	for m in ords do

		if m le TransitiveGroupDatabaseLimit() then

			First:=TransitiveGroups(m);
			Second:=[G : G in First | Order(G) eq m];
			Third:=[G: G in Second | IsSolvable(G)];

			T:=Third;
		else

			Sm:=Sym(m);

			T:=Subgroups(Sm:IsTransitive:=true, IsSolvable:=true, OrderEqual:=m);
			T:=[r`subgroup: r in T];
		end if;

		possi:=possi cat T;

	end for;

	return possi;

end function;



function SinclairStartingGroup(vertexlist,slopelist, residualpolynomials, BaseField,phi_0)

	vertices:=vertexlist;
	slopes:=slopelist;
	
	K:=BaseField;
	Ky<y>:=PolynomialRing(K);
	pi_K:=UniformizingElement(K);

	n:=RecoverPolyDegree(vertices);
	p:=Prime(PrimeRing(K));
	m:=Valuation(n,p);

	if #slopes eq 1 then
		return Ram_OneSegment(slopes[1], residualpolynomials[1],K,n,1);
	end if;


	//Find e0.

	e0:=n/(p^m);
	e0:=Integers()!e0;

	c:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices

	prod:=Sym(1);  //initial factor
	WreathFactors:=[];


	for i in [1..#slopes] do
		PolyDegree:=Integers()!(c[i+1]/c[i]);	//Degree is ratio of x-coord of consecutive endpoints.
		leftend:=c[i];

		if slopes[i] eq 0 then

			if e0 ne PolyDegree then
				error "Error: Tame degree is wrong.";
			end if;

			r,a,b:=Xgcd(e0,1);

			tau:=y^e0 - (K!(-phi_0/pi_K)^b) * pi_K;

			temp:=TameGaloisGroup(ext<K|tau>);  

		else
				
			temp:=Ram_OneSegment(slopes[i], residualpolynomials[i],K,PolyDegree,leftend);

		end if;

		Append(~WreathFactors,temp);
	end for;


	WreathFactors:=Reverse(WreathFactors);

	for i in [1..#WreathFactors] do
		Factt:=WreathFactors[i];
		prod:=WreathProduct(Factt,prod);
	end for;

	return prod;

end function;



function Ram_MaxTame(Vertexlist,slopelist,residual_polynomials,BaseField, phi_0, resexp)

	//resexp:  number that roots of residual polynomials must be raised to.

	vertices:=Vertexlist;
	slopes:=slopelist;

	if vertices[1][1] eq 0 then				//may not be needed.
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	K:=BaseField;	
	
	p:=Prime(PrimeRing(K));
	n:=RecoverPolyDegree(vertices);	

	m:=Valuation(n,p);
	e0:=n div p^m;
	e0:=Integers()!e0;

	if e0 eq 1 then
		l:=#slopes;
	else  //e0>1
		l:=#slopes-1;		//l+1 st segment is horizontal on x-axis
	end if;

	a:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices

	s:=[];

	for j in [1..(l+1)] do
		s[j]:=Valuation(a[j],p);
	end for;

	e:=[];
	h:=[];

	for i in [1..l] do
		e[i]:=Denominator(-slopes[i]);		//list of (negative) slope denominators
		h[i]:=Numerator(-slopes[i]);		//list of (negative) slope numerators
	end for;

	gamma:=[];
	
	f:=[];						//Segmental Inertia of segments

	RK,sigma:=ResidueClassField(K);
	Rs<r>:=PolynomialRing(RK);
	
	A:=residual_polynomials;
	A:=A[1..l];

	for i in [1..l] do

		f[i]:=Segmental_Inertia(A[i]);
		
	end for;

	Indicess:=[UnityIndex(K,e0*e[i]): i in [1..l]];
	F:=LCM(Flat([f,Indicess]));

	I:=UnramifiedExtension(K,F);    

	RI,tau:=ResidueClassField(I);
	Rw<w>:=PolynomialRing(RI);

	for i in [1..l] do

		Ai:=A[i];
		Ai:=Rs!Ai;			//Might be unnecessary, magma knows that RL=RK

		Ai:=Rw!Ai;
		gamma_i:=Roots(Ai:Max:=1);
		gamma_i:=gamma_i[1][1]^resexp;		//note gamma_i lives in RI
		gamma_i:=(gamma_i @@ tau);	//find pre_image of gamma_i, this is what we want
		Append(~gamma,gamma_i);

	end for;

	Iz<z>:=PolynomialRing(I);

	Tame_extensions:=[]; 

	for i in [1..l] do
		
		_, bi ,di:=Xgcd(h[i],e[i]);
		vi:=e0 * p^(m-s[i]) + n +1;		

		Ti:=TotallyRamifiedExtension(I, z^(e0*e[i]) - ((-1)^vi)* phi_0*gamma[i]^(bi*n) );
		Append(~Tame_extensions,Ti);
	end for;

	Towe:=TameCompositum(Tame_extensions);

	return Towe;	

end function;




function Sinclair_TwoSegmentCase(Gal1,N1, vertices, slopes,residualpolynomials,BaseField,phi_0, resexp:A:=0, PParity:="")
	//resexp:  the number that roots of residual polynomials must be raised to when computing T/K

	resexp:=Integers()!resexp;

	if #slopes ne #residualpolynomials then
		error "Error:  Must have a residual polynomial for each segment";  //Catch possible recursion error.
	end if;

	K:=BaseField;
	Ky<y>:=PolynomialRing(K);
	pi_K:=UniformizingElement(K);

	p:=Prime(PrimeRing(K));
	n:=RecoverPolyDegree(vertices);

	//Gal1 is Gal(Li-1/Li) for segment being included now.
	//N1 is Gal(Li/K)

	W:=WreathProduct(Gal1,N1);

	Num_Remaining:=[];

	if vertices[1][1] eq 1 then
		ps1:=vertices[2][1];
	else
		ps1:=vertices[1][1];
	end if;

	s1:=Valuation(ps1,p);
	T:=Ram_MaxTame(vertices,slopes,residualpolynomials,K, phi_0, resexp);

	//Look at information of rp(L1/K) to get T1, the max tame subext of L1/K.
	
	lvert, lslope, shrinkfactor:=NextLowerPolygon(vertices, slopes);
	shrinkfactor:=Integers()!shrinkfactor;
	shrinkfactor:=shrinkfactor*resexp;

	if lslope eq [0] then  //second segment is horizontal so T1 is L1
		m:=Valuation(n,p);
		e0:=n/(p^m);
		e0:=Integers()!e0;

		r,a,b:=Xgcd(e0,1);

		tau:=y^e0 - (K!(-phi_0/pi_K)^b) * pi_K;
		T1:=ext<K|tau>;
		GalT1:=N1;

	else
		T1:=Ram_MaxTame(lvert,lslope,residualpolynomials[2..#residualpolynomials], K, phi_0,shrinkfactor);
		GalT1:=TameGaloisGroup(T1);
	end if;

	//Find DL1 which is Degree(L1,K).

	DL1:=RecoverPolyDegree(lvert);

	if DL1 ne Integers()!(n/ps1) then
		error "Error: DL1 is incorrect.";
	end if;

	//Start filtering.

	Pos_Order:=PossibleCandidateOrders(N1,s1,T,T1,DL1);

	if n lt 35 then
vprint Galois,1:"GaloisGroupMilstead: Finding transitive subgroups of", W;

		Can:=Subgroups(W:IsTransitive:=true);
		C:=[r`subgroup: r in Can];

		First:=Keys(CandidateClasses(C));
		Num_Remaining[1]:=#First;
vprint Galois,1:"GaloisGroupMilstead: After 1", Num_Remaining;

		Second:=[c: c in First| Order(c) in Pos_Order];
		Num_Remaining[2]:=#Second;
	else

		min_order:=Min(Pos_Order);
		max_order:=Max(Pos_Order); //Use list of possible orders to set parameters for calculating initial group candidates as subgroups of W.

		max_order:=Min(max_order,Order(W));  //Need to limit ourselves to size of W.
vprint Galois,1:"GaloisGroupMilstead: Finding transitive subgroups wof", W,"whose order divides",max_order;
		Can:=Subgroups(W:IsTransitive:=true, OrderMultipleOf:=min_order, OrderDividing:=max_order );
		C:=[r`subgroup: r in Can];
		InitialCandidates:=CandidateClasses(C);  //associative array

		First:=Keys(InitialCandidates);
		Num_Remaining[1]:=#First;
vprint Galois,1:"GaloisGroupMilstead: After 1", Num_Remaining;
		Second:=First;
		Num_Remaining[2]:=#Second;

	end if;

vprint Galois,1:"GaloisGroupMilstead: After 2", Num_Remaining;

	if Type(PParity) eq MonStgElt then   //PParity not specified.

		ValD:=RecoverDiscriminant(vertices);
		ValD:=Integers()!ValD;

		if IsOdd(ValD) then
			PParity:={-1};
		else
			PParity:={1,-1};
		end if;

		third:=[G: G in Second | (GroupParity(G) in PParity)];

	elif PParity in {1, -1} then

		third:=[G: G in Second | (PParity eq GroupParity(G))];
	else

		error "Error: PParity is of the wrong type.";
	end if;


	Num_Remaining[3]:=#third;

vprint Galois,1:"GaloisGroupMilstead: After 3", Num_Remaining;
	if #third le 1 then
		return third, Num_Remaining;
	end if;

	Fourth:=[];

	if Type(A) eq RngIntElt then  //Automorphism group not specified

		Autoo:=Possible_Auto_Orders(vertices, slopes, residualpolynomials);

		for i in [1..#third] do
			G:=third[i];
			Cen:=Centralizer(Sym(n),G);

			if Order(Cen) in Autoo then
				Append(~Fourth,G);
			end if;

		end for;

	else

		Autoo:=A;

		for i in [1..#third] do
			G:=third[i];
			Cen:=Centralizer(Sym(n),G);
	
			if IsIsomorphic(Autoo,Cen) then
				Append(~Fourth,G);
			end if;

		end for;

	end if;

	Num_Remaining[4]:=#Fourth;

vprint Galois,1:"GaloisGroupMilstead: After 4:", Num_Remaining;

	if #Fourth le 1 then
		return Fourth, Num_Remaining;
	end if;

	BTuples:=[];

	for n in Fourth do
		bool,B1:=PossibleB1(n,GalT1);
		
		if bool then
			Append(~BTuples,<n,B1>);
		end if;
	end for;

	TameGrp:=TameGaloisGroup(T);
	Fifth:=[];  //after using in Num_Remaining, rename it BSTuples

	for m in BTuples do  //now over BTuples
		
		Candidate:=m[1];

		B:=m[2];  		//Possible B's.  Loop over this.

		BS_set:={};  //Set of tuples <B1,[S]>
		for B1 in B do

			bool,S:=AMP_Max_Tame_Test(Candidate, TameGrp, p, B1);
			//bool;
	
			if bool then
				BS_set:=BS_set join {<B1,S>};
			end if;
		end for;
		
		if #BS_set gt 0 then
	
			Append(~Fifth,<Candidate,BS_set>);
		end if;

	end for;

	Num_Remaining[5]:=#Fifth;

vprint Galois,1:"GaloisGroupMilstead: After 5:", Num_Remaining;
	
	if #Fifth le 1 then

		return [t[1]: t in Fifth], Num_Remaining;
	end if;


	BSTuples:=Fifth;

	//return Num_Remaining;

	Sixth:=[];  	//after using in Num_Remaining, rename it BSU_Tuples.  Contains B1, S, and U.

	for m in BSTuples do
		Candidate:=m[1];
		CurrentSet:=m[2];  //of the form {<B1, [S]>}.

vprint Galois,1:"GaloisGroupMilstead: #B1",#CurrentSet;

		BSU_set:={};

		for C in CurrentSet do
			B1:=C[1];
			Possible_S:=C[2];

			bool,U:=PrevNormalTest(Candidate,N1,B1);

			if bool then
				BSU_set:= BSU_set join {<B1,Possible_S,U> };
			end if;	

		end for;

		if #BSU_set gt 0 then
			Append(~Sixth,<Candidate,BSU_set>);
		end if;

	end for;

	Num_Remaining[6]:=#Sixth;

vprint Galois,1:"GaloisGroupMilstead: After 6:", Num_Remaining;

	if #Sixth le 1 then
		return [s[1]: s in Sixth], Num_Remaining;
	end if;

	BSU_Tuples:=Sixth;

	//Change the next four lines

	Seventh:=[];

	for t in BSU_Tuples do

		H:=t[1];
		CurrentSet:=t[2];

		New_set:={};   // {<B1, SU_set>}

		for C in CurrentSet do
			B1:=C[1];
			Possible_S:=C[2];
			Possible_U:=C[3];
			SU_set:={};			//{<s,u>}

			for s in Possible_S do

				for u in Possible_U do
					bool:=IsIsomorphic(quo<B1|s>, quo<u|(s meet u)>);
				

					if bool then
						SU_set:=SU_set join {<s,u>};
					end if;

				end for;

			end for;

			if #SU_set gt 0 then
				New_set:=New_set join {<B1, SU_set>};

			end if;


		end for;

		if #New_set gt 0 then
			Append(~Seventh,<H,New_set>);
		end if;

	end for;

	Num_Remaining[7]:=#Seventh;

vprint Galois,1:"GaloisGroupMilstead: After 7:", Num_Remaining;

	if #Seventh le 1 then
		return [s[1]: s in Seventh], Num_Remaining;
	end if;
	
	Eighth:=[];

	for v in Seventh do

		H:=v[1];
		Btupleset:=v[2];

		bbool, D:=PossibleD(H,DL1,Gal1);

		if bbool then
			Append(~Eighth,<H,Btupleset,D>);
		end if;

	end for;
	

	Num_Remaining[8]:=#Eighth;

vprint Galois,1:"GaloisGroupMilstead: After 8:", Num_Remaining;

	if #[r[1]: r in Eighth] le 1 then
		return [r[1]: r in Eighth], Num_Remaining;
	end if;
	
	Ninth:=[];

	for r in Eighth do

		H:=r[1];

		BB:={};   //will be {<B1,s,u,D1,D0>}

		for T in r[2] do
			B1:=T[1];
			su_pairs:=T[2];

			for w in su_pairs do
				s:=w[1];
				u:=w[2];

				for z in r[3] do
		        		D1:=z[1];
                        		if u subset D1 then
						for D0 in z[2] do
							Group1:=quo<(s meet D1)|(s meet D0)>;
							bool1:= (Order(Group1) eq p^s1) and IsElementaryAbelian(Group1);
							Group2:=quo<(s meet u)|((s meet u) meet D0)>;
							bool2:=(Order(Group2) le p^s1) and IsElementaryAbelian(Group2) and IspGroup(Group2,p);

							TN1L0_TN1 := #(u meet s)/#(D0 meet u meet s);
							TN1_K := #H/#(u meet s);
							e0 := DL1/p^Valuation(DL1,p);
							N1_T1 := Order(N1)/Degree(T1,K);
							possible_orders := {TN1_K*p^w:w in [Valuation(TN1L0_TN1,p)..e0*N1_T1*Valuation(TN1L0_TN1,p)]};

                                        		bool3:=#H in possible_orders;
										
							if (bool1 and bool2 and bool3) then
								BB:=BB join {<B1,s,u,D1,D0>};
							end if;

			 			end for;

                            		end if;
				end for;

			end for;
		end for;

		if #BB gt 0 then
			Append(~Ninth,<H,BB>);
		end if;
	end for;

	Num_Remaining[9]:=#Ninth;
vprint Galois,1:"GaloisGroupMilstead: After 9:", Num_Remaining;
	 
	Candidates:=[Ninth[i][1] : i in [1..#Ninth]];

	return Candidates, Num_Remaining;

end function;



function Sinclair_ThreeSegmentCase(vertices, slopes,residualpolynomials,BaseField,phi_0, WreathFactors:A:=0,PParity:="")


	//We first do segments 2 and 3.  The first is incorporated after.	

	lvert, lslope, shrinkfactor:=NextLowerPolygon(vertices, slopes);

	PossibleN1,_:=Sinclair_TwoSegmentCase(WreathFactors[2],WreathFactors[3],lvert, lslope, residualpolynomials[2..#residualpolynomials], BaseField,phi_0, shrinkfactor);


	//Now incorporate segment 1

	Gal1:=WreathFactors[1];
	Candidates:=[];

vprint Galois,1:"GaloisGroupMilstead: Num N1", #PossibleN1;

	if Type(PParity) eq MonStgElt then   //PParity not specified.

		for N1 in PossibleN1 do

vprint Galois,1:"GaloisGroupMilstead: N1", TransitiveGroupIdentification(N1);

			if Type(A) eq RngIntElt then  //Automorphism group not specified

				a,b:=Sinclair_TwoSegmentCase(Gal1,N1, vertices, slopes,residualpolynomials,BaseField,phi_0, 1); 
				Candidates:=Candidates cat a;
			else

				Autoo:=A;

				a,b:=Sinclair_TwoSegmentCase(Gal1,N1, vertices, slopes,residualpolynomials,BaseField,phi_0, 1:A:=Autoo);
				Candidates:=Candidates cat a;

			end if;

		end for;

	elif PParity in {1,-1} then
		PP:=PParity;

		for N1 in PossibleN1 do

vprint Galois,1:"GaloisGroupMilstead: N1", TransitiveGroupIdentification(N1);

			if Type(A) eq RngIntElt then  //Automorphism group not specified

				a,b:=Sinclair_TwoSegmentCase(Gal1,N1, vertices, slopes,residualpolynomials,BaseField,phi_0, 1:PParity:=PP); 
				Candidates:=Candidates cat a;
			else

				Autoo:=A;

				a,b:=Sinclair_TwoSegmentCase(Gal1,N1, vertices, slopes,residualpolynomials,BaseField,phi_0, 1:A:=Autoo, PParity:=PP);
				Candidates:=Candidates cat a;

			end if;

		end for;

	else 
		error "Error: PParity is of the wrong type.";
	end if;	


	if IsEmpty(Candidates) then 	//Can happen with bad automorphism group

		return [];
	end if;

	//sort by trans id.

	return SetToSequence(Keys(CandidateClasses(Candidates)));

end function; 
	

function SinclairMain(vertexlist, residualpolynomials, BaseField,phi_0:A:=0,PParity:="")
//Primary function

	//1. Process inputs

	vertices:=ProcessVertices(vertexlist);	

	slopes:=[];

	for i in [1..(#vertices -1)] do  	//Find slopes using the function MySlope

		tup1:=vertices[i];
		tup2:=vertices[i+1];
		s:=MySlope(tup1,tup2);

		Append(~slopes,s);
	end for;

	K:=BaseField;
	Ky<y>:=PolynomialRing(K);
	pi_K:=UniformizingElement(K);

	n:=RecoverPolyDegree(vertices);

	p:=Prime(PrimeRing(K));
	m:=Valuation(n,p);

	e0:=n/(p^m);
	e0:=Integers()!e0;


	if Type(residualpolynomials) eq Tup then

		residualpolynomials:=TupleToSequence(residualpolynomials);
	end if;


	//2. Early exits for one segment case.

	if n eq e0 then

		r,a,b:=Xgcd(e0,1);
		tau:=y^e0 - (K!(-phi_0/pi_K)^b) * pi_K;

		return [TameGaloisGroup(ext<K|tau>)];

	end if;

	
	if #slopes eq 1 then
		G:=Ram_OneSegment(slopes[1], residualpolynomials[1],K,n,1);

		return [G];
	end if;


	//3. Find Gal(L_{i-1}/L_i) for each i store in list.

	c:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices

	WreathFactors:=[];

	for i in [1..#slopes] do
		PolyDegree:=Integers()!(c[i+1]/c[i]);	//Degree is ratio of x-coord of consecutive endpoints.
		leftend:=c[i];

		if slopes[i] eq 0 then

			r,a,b:=Xgcd(e0,1);
			tau:=y^e0 - (K!(-phi_0/pi_K)^b) * pi_K;

			temp:=TameGaloisGroup(ext<K|tau>);  

		else
				
			temp:=Ram_OneSegment(slopes[i], residualpolynomials[i],K,PolyDegree,leftend);

		end if;

		Append(~WreathFactors,temp);
	end for; 
	//Currently WreathFactors is [Gal(L0/L1), Gal(L1/L2), Gal(L2/L3), ..., Gal(Ll/K)]


	//4. Two Segment Case

	if #slopes eq 2 then


		if Type(PParity) eq MonStgElt then	//PParity not specified.


			if Type(A) eq RngIntElt then  //Automorphism group not specified

				return Sinclair_TwoSegmentCase(WreathFactors[1],WreathFactors[2],vertices, slopes,residualpolynomials,K,phi_0,1);
			else

				Autoo:=A;

				return Sinclair_TwoSegmentCase(WreathFactors[1],WreathFactors[2],vertices, slopes,residualpolynomials,K,phi_0,1:A:=Autoo);
			end if;


		elif PParity in {1,-1} then

			PP:=PParity;

			if Type(A) eq RngIntElt then  //Automorphism group not specified

				return Sinclair_TwoSegmentCase(WreathFactors[1],WreathFactors[2],vertices, slopes,residualpolynomials,K,phi_0,1:PParity:=PP);
			else

				Autoo:=A;

				return Sinclair_TwoSegmentCase(WreathFactors[1],WreathFactors[2],vertices, slopes,residualpolynomials,K,phi_0,1:A:=Autoo, PParity:=PP);
			end if;
		else
			error "Error:  PParity is of the wrong type";
		end if;
		 
	end if;


	if #slopes eq 3 then 


		if Type(PParity) eq MonStgElt then	//PParity not specified.

			if Type(A) eq RngIntElt then  //Automorphism group not specified


				return 	Sinclair_ThreeSegmentCase(vertices, slopes,residualpolynomials,BaseField,phi_0, WreathFactors);
			else

				Autoo:=A;

				return  Sinclair_ThreeSegmentCase(vertices, slopes,residualpolynomials,BaseField,phi_0, WreathFactors:A:=Autoo);
			end if;


		elif PParity in {1,-1} then

			PP:=PParity;

			if Type(A) eq RngIntElt then  //Automorphism group not specified


				return 	Sinclair_ThreeSegmentCase(vertices, slopes,residualpolynomials,BaseField,phi_0, WreathFactors:PParity:=PP);
			else

				Autoo:=A;

				return  Sinclair_ThreeSegmentCase(vertices, slopes,residualpolynomials,BaseField,phi_0, WreathFactors:A:=Autoo, PParity:=PP);
			end if;

		else
			error "Error:  PParity is of the wrong type";
		end if;
		 
	end if;


	if #slopes gt 3 then
		error "Error: This case isn't implemented yet.";
	end if;

	return 0;

end function;





























//Verification functions below.


function My_ResidualPolynomials(ramification_polygon,rho)

	//rho is ramification polynomial.

	L:= CoefficientRing(rho);

	vertices, slopes:=Recover_Vertices_Slopes(ramification_polygon);

	if vertices[1][1] eq 0 then
		vertices:=vertices[2..#vertices];		//get rid of point with x-coordinate=0
		slopes:=slopes[2..#slopes];			//get rid of segment with infinite slope
	end if;

	a:=[Integers()!vertices[i][1]: i in [1..#vertices]];	//list of x-coordinates of vertices
	b:=[Integers()!vertices[i][2]: i in [1..#vertices]];	//list of y-coordinates of vertices

	l:=#slopes;

	pi_L:=UniformizingElement(L);
	
	e:=[];
	h:=[];

	for i in [1..l] do
		e[i]:=Denominator(-slopes[i]);		//list of (negative) slope denominators
		h[i]:=Numerator(-slopes[i]);			//list of (negative) slope numerators
	end for;

	A:=[];
	
	RL,omega:=ResidueClassField(L);
	Ry<y>:=PolynomialRing(RL);

	for i in [1..l] do
		di:=a[i+1]-a[i];
		fin:=Integers()!(di/e[i]);
		Ai:=Polynomial(Ry,[&+[omega(Coefficient(rho,j*e[i]+a[i]) div pi_L^(-j*h[i]+b[i]))*y^j: j in [0..fin]]]);
		Ai:=Ry!Ai;
		A[i]:=Ai;
	end for;


	return A;


end function;



function Check_LowerPolygon(phi)

	if Number_Of_Segments(phi) eq 1 then
		return true;
	end if;
	
	K:=CoefficientRing(phi);
	stem:=ext<K|phi>;

	L:=Ramification_Segments(stem);
	 L:=L[#L];

	//First thing for comparison is the output of NextLowerPolygon

	ramification_polygon1,rho1:=My_RamificationPolygon(phi);
	vertices1, slopes1:=Recover_Vertices_Slopes(ramification_polygon1);

	Firstv, Firsts:=NextLowerPolygon(vertices1, slopes1);

	//Second thing for comparison is rampolygon of L1/K.

	L1:=BaseRing(L);

	psi:=DefiningPolynomial(L1,PrimeRing(L1));

	ramification_polygon2, rho2:=My_RamificationPolygon(psi);
	vertices2, slopes2:=Recover_Vertices_Slopes(ramification_polygon2);

	//Comparisons

	bool1:=Firstv eq vertices2;
	bool2:=Firsts eq slopes2;

	bool3:=Valuation(Discriminant(psi)) eq RecoverDiscriminant(vertices2);


	if {bool1, bool2, bool3} eq {true} then

		return true;
	else
		return false;
	end if;

	return false;

end function;


function VerifyMain(phi,Galid)

	//Galid is the transitive group Identification of Gal(phi).

	if IsEisenstein(phi) eq false then
		error "Error: inputted polynomial isn't Eisenstein";
	end if;

	ramification_polygon,rho:=My_RamificationPolygon(phi);
	vertices, slopes:=Recover_Vertices_Slopes(ramification_polygon);

	prec:=Precision(CoefficientRing(phi));

	my_discriminant:=RecoverDiscriminant(vertices);        //Find valuation of discriminant.
	my_discriminant:=Integers()!my_discriminant;

	if prec le my_discriminant then

    		error "Error: Precision of polynomial is too low.";
       
	end if;

	residualpolynomials:=My_ResidualPolynomials(ramification_polygon,rho);

	BaseField:=CoefficientRing(phi);

	phi_0:=ConstantCoefficient(phi);

	if #slopes eq 1 then

		Candidates:=SinclairMain(vertices, residualpolynomials, BaseField,phi_0);

		Cand_Num:=[TransitiveGroupIdentification(r): r in Candidates];

		return Galid in Cand_Num;

	end if;		

	Autoo:=AutomorphismGroup(ext<BaseField|phi>,BaseField);	//Currently we are specifying automorphism group.
	
	Candidates:=SinclairMain(vertices, residualpolynomials, BaseField,phi_0:A:=Autoo, PParity:=PolynomialParity(phi));

	Cand_Num:=[TransitiveGroupIdentification(r): r in Candidates];

	return Galid in Cand_Num;

end function;


function FilteringPower(phi)

	
	if IsEisenstein(phi) eq false then
		error "Error: inputted polynomial isn't Eisenstein";
	end if;

	ramification_polygon,rho:=My_RamificationPolygon(phi);
	vertices, slopes:=Recover_Vertices_Slopes(ramification_polygon);

	if #slopes gt 3 then
		error "Error: Not implemented yet";
	end if;

	prec:=Precision(CoefficientRing(phi));

	my_discriminant:=RecoverDiscriminant(vertices);        //Find valuation of discriminant.
	my_discriminant:=Integers()!my_discriminant;

	if prec le my_discriminant then

    		error "Error: Precision of polynomial is too low.";
       
	end if;


	residualpolynomials:=My_ResidualPolynomials(ramification_polygon,rho);

	BaseField:=CoefficientRing(phi);

	phi_0:=ConstantCoefficient(phi);

	if #slopes eq 1 then

		return SinclairMain(vertices, residualpolynomials, BaseField,phi_0);
	end if;
	
	Autoo:=AutomorphismGroup(ext<BaseField|phi>,BaseField);  //Currently we are specifying automorphism group.
	
	//Candidates, Num_Remaining:=SinclairMain(vertices, residualpolynomials, BaseField,phi_0:A:=Autoo, PParity:=PolynomialParity(phi));

	//return Num_Remaining;

	return SinclairMain(vertices, residualpolynomials, BaseField,phi_0:A:=Autoo, PParity:=PolynomialParity(phi));


end function;



function ResidualTheory(phi,order)

//Currently scrapped.  Here is a brief description of the idea.

//If the degree of phi is a power of a prime p and all the residual polynomials are square free then
// Gal(phi) has order Degree(phi)* [T:K].

	if IsEisenstein(phi) eq false then
		error "Error: inputted polynomial isn't Eisenstein";
	end if;

	ramification_polygon,rho:=My_RamificationPolygon(phi);
	vertices, slopes:=Recover_Vertices_Slopes(ramification_polygon);

	residualpolynomials:=My_ResidualPolynomials(ramification_polygon,rho);

	Squarebool:=[IsitSquareFree(r): r in residualpolynomials];

	if Set(Squarebool) eq {true} then

		K:=CoefficientRing(phi);
		T:=Max_Tame_Subextension(phi);

		Comparr:=Degree(phi)*Degree(T,K);

		return Comparr eq order;
	end if;

	return false;

end function;	




function DetermineMySlopes(vertexlist)

//Given the list of vertices of a ramification polygon, returns the list of slopes of the polygon.

	vertices:=ProcessVertices(vertexlist);	

	slopes:=[];

	for i in [1..(#vertices -1)] do  	//Find slopes using the function MySlope

		tup1:=vertices[i];
		tup2:=vertices[i+1];
		s:=MySlope(tup1,tup2);

		Append(~slopes,s);
	end for;

	return slopes;

end function;

function CompleteTime(t)


	dayq:=t div 86400;
	dayr:=t -86400*dayq;

	hourq:= dayr div 3600;
	hourr:=dayr - 3600*hourq;

	minq:=hourr div 60;
	minr:=hourr-60*minq;

	return dayq, "Days", hourq, "Hours", minq, "Minutes", minr, "Seconds";

end function;



function MyPoliFact(phi)

	//Return the exact factors of phi over its stem field L.	

	K:=CoefficientRing(phi);
	L:=ext<K|phi>; //StemField
	Ly<y>:=PolynomialRing(L);

	phi:=Ly!phi;

	FF:=Factorization(phi);
	F:=[f[1]: f in FF];		//Don't need the multiplicities.

	return F;

end function;

intrinsic GaloisGroupMilstead(phi:global:=false) -> .
{}
   if not IsEisenstein(phi) then
     fact, c, exts := Factorization(phi:Extensions:=true);
     if #fact gt 1 then
        error "GaloisGroupMilstead: Polynomial must be irreducible";
     end if;
     if exts[1]`F gt 1 then
        error "GaloisGroupMilstead: Polynomial must generate totally ramified extension";
     end if;
     phi := DefiningPolynomial(exts[1]`Extension);
   end if;

   return Eisenstein_Case(phi:global:=global);

end intrinsic;

