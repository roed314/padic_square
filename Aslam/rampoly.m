//RamificationPolygon: f(x) eisenstein -> Knoten des Verzweigungspolygons von f(x)
//Flatten: f(x) -> assoziierte Polynome zu den Segmenten des Newton-Polygons von f(x)
//FlatRamificationPolygon: f(x) eisenstein -> zu jedem Segment des Verzweigungspolygons: < Nenner der Steigung, ass. Polynom >
// + verschiedene Varianten der obigen Verfahren

///////////////////////////////////////////////////////////////////////////////////////////
//vgl. Definition 4.6
intrinsic RamificationPolynom(f::RngUPolElt) -> RngUPolElt
	{f(x) eisenstein -> Verzweigungspolynom vom f(x)}

	Qp:=CoefficientRing(f);
	L<alpha>:=TotallyRamifiedExtension(Qp,f);
	P<x>:=PolynomialRing(L);
	rho:=Evaluate(f,alpha*x+alpha) div (x*alpha^Degree(f));
       
        return rho,L;

end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////
// Berechnet das Verzweigungspolygon eines Eisensteinpolynoms mit Korollar 4.11

intrinsic RamificationPolygon(f::RngUPolElt) -> SeqEnum
{Computes the RamificationPolygon of an Eisenstein polynomial}

	require LeadingCoefficient(f) eq 1: "Polynomial is not monic.";  
	assert2 IsEisenstein(f); 

	n := Degree(f);
	K := CoefficientRing(f);
	p := Prime(K);
	e_K := RamificationDegree(K);
	r := Valuation(n,p);

	Val := [n*Valuation(Coefficient(f,j)) : j in [0..n] ];
	L := [];

	for i in [0..r-1] do	// valuation of p^i coefficients
		valbinom := 0;
pts := [* (Val[p^i+1] + p^i - n) *]; //16 Nov 2017 we change it because f:=x^9+9*x+3; didn't work.
		for j in [p^i+1..n] do // recursionformula to compute valuation of binomial coefficients
			valbinom := Valuation(j,p) - Valuation(j - p^i,p) + valbinom;
			Append(~pts, e_K*n*valbinom + Val[j+1] + j - n);
		end for;
Append(~L,<p^i-1,Min([pts[i] : i in [1..#pts]])>);//16 Nov 2017 we change it because f:=x^9+9*x+3; didn't work.
	end for;

	Append(~L,<p^r-1,0>); // valuation always zero 
	k := 1;
	V := [L[1]];

	while L[k][2] ne 0 do // exit condition needs monic poly
		delta := [(L[i][2]-L[k][2])/(L[i][1]-L[k][1]) : i in [k+1..#L]];
		min := Min(delta);
		j := Index(delta, min) + k;
		h := Index(delta, min, j-k+1);

		while h ne 0 do
			j := h+k;
			h := Index(delta, min, j-k+1);
		end while;

		Append(~V, L[j]);
		k := j;
	end while;

	if n-1 ne V[#V][1] then
		Append(~V,<n - 1,0>);
	end if;

	return V;
end intrinsic; 

///////////////////////////////////////////////////////////////////////////////////////
//vgl. Definition 4.1
intrinsic Flatten(f::RngUPolElt) -> SeqEnum
	{Berechnung der assoziierten Polynome zu den Segmenten des Newton-Polygons von f(x).}

	v:=Vertices(NewtonPolygon(f));
	Polys:=[];

	for i in [2..#v] do
		E:=Integers() ! (v[i][1]-v[i-1][1]);
		H:=Integers() ! (v[i-1][2]-v[i][2]);
		d:=GCD(E,H);
		e:=E div d;
		h:=H div d;
		K:=CoefficientRing(f);
		p:=UniformizingElement(K);
		R:=ResidueClassField(Integers(K));
		RX<x>:=PolynomialRing(R);
		gi:=Zero(RX);
		
		for j in [0..d] do
			a:=Integers() ! -v[i-1][2];
			b:=Integers() ! v[i-1][1];
			aj:=R ! (p^(a+j*h)*Coefficient(f,b+j*e));
			gi:=gi+aj*x^j;
		end for;
		
		Append(~Polys,<e,gi>);
	end for;

	return Polys;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////////////////////////
//wie Flatten, es werden aber mehr Informationen ausgegeben (fuer pReduction !!)
intrinsic XFlatten(f::RngUPolElt) -> SeqEnum
	{}

	v:=Vertices(NewtonPolygon(f));
	Polys:=[];

	for i in [2..#v] do
		E:=Integers() ! (v[i][1]-v[i-1][1]);
		H:=Integers() ! (v[i-1][2]-v[i][2]);
		d:=GCD(E,H);
		e:=E div d;
		h:=H div d;
		K:=CoefficientRing(f);
		p:=UniformizingElement(K);
		R:=ResidueClassField(Integers(K));
		RX<x>:=PolynomialRing(R);
		gi:=Zero(RX);
		
		for j in [0..d] do
			a:=Integers() ! -v[i-1][2];
			b:=Integers() ! v[i-1][1];
			aj:=R ! (p^(a+j*h)*Coefficient(f,b+j*e));
			gi:=gi+aj*x^j;
		end for;
		
		Append(~Polys,< h, e, (Integers()!v[i-1][1]+1), gi >);	//ÄNDERUNG "h" und "v[i-1][1]+1"=p^{s-{i-1}}
	end for;

	return Polys;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////////////
//wie Flatten, aber das ass. Polynom über K, nicht über RK
intrinsic Flatten2(f::RngUPolElt) -> SeqEnum  
	{}

	v:=Vertices(NewtonPolygon(f));
	Polys:=[];

	for i in [2..#v] do
		E:=Integers() ! (v[i][1]-v[i-1][1]);
		H:=Integers() ! (v[i-1][2]-v[i][2]);
		d:=GCD(E,H);
		e:=E div d;
		h:=H div d;
		K:=CoefficientRing(f);
		p:=UniformizingElement(K);
		R:=ResidueClassField(Integers(K));
		KX<x>:=PolynomialRing(K);	
		//RX<x>:=PolynomialRing(R);
		//gi:=Zero(RX);
		gi:=Zero(KX);
		
		for j in [0..d] do
			a:=Integers() ! -v[i-1][2];
			b:=Integers() ! v[i-1][1];
			aj:=(p^(a+j*h)*Coefficient(f,b+j*e));   //liefert polynom über K!
			//aj:=R ! (p^(a+j*h)*Coefficient(f,b+j*e));
			gi:=gi+aj*x^j;
		end for;
		
		Append(~Polys,<e,gi>);
	end for;

	return Polys;

end intrinsic;

////////////////////////////////////////////////////////////////////////////////////////////
// Compute associated polynomials belonging to the parts of f's Ramificationpolygon.
// xflat gives more information, required for pReduction 
intrinsic FlatRamificationPolygon(f::RngUPolElt:xflat:=false) -> SeqEnum
  {f(x) eisenstein -> zu jedem Segment des Verzweigungspolygons: < Nenner der Steigung, ass. Polynom >}
 
  require LeadingCoefficient(f) eq 1: "Polynomial is not monic.";
  assert2 Type(xflat) eq BoolElt;
  assert2 IsEisenstein(f); 

  n := Degree(f);
  K := CoefficientRing(f);
  p := Prime(K);
  r := Valuation(n,p);
  e_K := RamificationDegree(K);

  Val := [n*Valuation(Coefficient(f,j)) : j in [0..n] ];
  L := [];
  minIndex := [];
//  pts := [];
  for i in [0..r-1] do	// valuation of p^i coefficients
	valbinom := 0;
	pts :=[* *];
	Append(~pts,  Val[p^i+1] + p^i - n ); 
	for j in [p^i+1..n] do // recursionformula to compute valuation of binomial coefficients
		valbinom := Valuation(j,p) - Valuation(j - p^i,p) + valbinom;
		Append(~pts, e_K*n*valbinom + Val[j+1] + j - n);
	end for;
	pts := [x: x in pts];
	Append(~L,<p^i-1,Min(pts)>);
        Append(~minIndex, Index(pts, Min(pts)) + p^i-1);	       
  end for;

  Append(~L,<p^r-1,0>); // valuation always zero 
  Append(~minIndex, n);

  k := 1;
  V := [L[1]];

  while L[k][2] ne 0 do // exit condition needs monic poly
	delta := [(L[i][2]-L[k][2])/(L[i][1]-L[k][1]) : i in [k+1..#L]];
	min := Min(delta);
	j := Index(delta, min) + k;
	h := Index(delta, min, j-k+1);

	while h ne 0 do
		j := h+k;
		h := Index(delta, min, j-k+1);
	end while;

	Append(~V, L[j]);
	k := j;
  end while;

  if n-1 ne V[#V][1] then 
	Append(~V,<n - 1,0>);
	Append(~L,<n-1,0>);
	Append(~minIndex,n);
  end if;

  R_K := ResidueClassField(Integers(K));
  pi := UniformizingElement(K);
  delta0 := R_K ! (-pi div Coefficient(f,0));
  delta1 := R_K ! (p div pi^e_K);
  R<x> := PolynomialRing(R_K);  
  Polys := [];

  for k in [2..#V] do // computes the associated polies using 4.1 

   	u := V[k-1][1];
   	v := V[k-1][2];	
   	E := V[k][1]-V[k-1][1];
   	H := V[k-1][2]-V[k][2];
   	d := GCD(E,H);
   	e := E div d;
   	h := H div d;
   	g := Zero(R);
	
         for i in [0..d] do
	        x:= <u + i * e, v - i *h >;
		if x notin L then // point not on line 
   		       continue;
   	        end if; 
  
		j := minIndex[Index(L,x)]; // details @ 4.1
   		r_2 := Valuation(Coefficient(f,j));
   		b := Binomial(j, u+i*e+1);
   		r_3 := Valuation(b,p);
   		delta2 := R_K ! (Coefficient(f,j) div pi^r_2);
   		delta3 := R_K ! (b div p^r_3);
   		gi:= delta0^(r_3*e_K+r_2)*delta1^(r_3)*delta2*delta3;
   		g := g + gi* R.1^i;
   	end for;
	if xflat then // more info for pReduction
		Append(~Polys, <h,e,V[k-1][1]+1,g>);
	else 
		Append(~Polys, <e,g>);
	end if;
  end for;

return Polys; 
end intrinsic;



///////////////////////////////////////////////////////////////////////////////////////
//wie FlatRamificationPolygon, aber es werden mehr Informationen ausgegeben (für pReduction !!)
intrinsic XFlatRamificationPolygon(f::RngUPolElt) -> SeqEnum
{}
	return FlatRamificationPolygon(f:xflat:=true);
end intrinsic;



//FlatRamificationPolygon mit ass. Polynomen ueber K 
intrinsic FlatRamificationPolygon2(f::RngUPolElt) -> SeqEnum
	{}

	Qp:=CoefficientRing(f);
	L<alpha>:=TotallyRamifiedExtension(Qp,f);
	P<x>:=PolynomialRing(L);
	rho:=Evaluate(f,alpha*x+alpha) div (x*alpha^Degree(f));

	return Flatten2(rho);

end intrinsic;

intrinsic Eisenstein(f::RngUPolElt) -> RngUPolElt
	{}
	_,_,ext:=Factorization(f:Extensions:=true);
	
	return DefiningPolynomial(ext[1]`Extension);
end intrinsic;
