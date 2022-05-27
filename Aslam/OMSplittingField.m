declare verbose OM, 6;


intrinsic AbsoluteNorm(g::RngUPolElt) -> RngUPolElt
{The absolute norm of g}

  B := CoefficientRing(Parent(g));
  while B cmpne BaseRing(B) do
    g := Norm(g);
    B := BaseRing(B);
  end while;
  return g;
end intrinsic;


////////////////////////////////////
// elt
////////////////////////////////////


function phi_expansion(f,phi:length:=0)
  expansion := [];
  l := 0;
  repeat
    l +:= 1;
    a := f mod phi;
    Append(~expansion,a);
    f := (f-a) div phi;
  until f eq 0 or l eq length;
  return expansion;
end function;

intrinsic elt_recformat() -> .
{}
// either unit and exponent or coeffs and exponents is assigned
  return recformat
        < frame,
          unit,      // a unit in OK elt = unit * pi^valuation
          coeffs,
          exponents,
          valuation   // the valuation of the element
        >;
end intrinsic;

intrinsic elt_print(A)
{}
  if A`frame`is_first then
    printf "(%o)*%o^%o",A`unit,UniformizingElement(CoefficientRing(Parent(A`frame`Phi))),Integers()!A`valuation;
  else
    printf " [ ";
    for i in [1..#A`coeffs] do
      printf "(";
      elt_print(A`coeffs[i]);
      printf ")*";
      printf "phi_%o^%o",A`frame`depth,A`exponents[i];
      if i ne #A`coeffs then
        printf " * ";
      end if;
    end for;
    printf " ] ";
  end if;
end intrinsic;


intrinsic elt_with_val(frame,v: extval := false) -> .
{}
// return an elt with valuation val
// if extval is true the exponent of pi is returned as a separate value
  OK := CoefficientRing(Parent(frame`Phi));
//"OK",OK;
//UniformizingElement(OK);
//Parent(v);
  if frame`is_first then
    if extval then
      return rec< elt_recformat() | frame := frame, unit := One(OK), valuation := 0 >, Integers()!v;
    else
      return rec< elt_recformat() | frame := frame, unit := One(OK), valuation := Integers()!v >,0;
    end if;
  else
    d := frame`segment`E/frame`segment`Eplus;
    tmpval := v * d;
    tmpool := frame`segment`val * d;
    mud := Denominator(tmpool);
    rcr := ResidueClassRing(Integers()!mud);
    tmpval := tmpval * mud;
    tmpool := tmpool * mud;
    sol := Integers()!((rcr!tmpval)/(rcr!tmpool));
//"sol",sol;
    rem_val := v - sol * frame`segment`val;
    coeff, val := elt_with_val(prev_frame(frame),rem_val: extval := extval);
    return rec< elt_recformat() | frame := frame, exponents := [sol], coeffs := [coeff] >, val;
  end if;
end intrinsic; 

intrinsic elt_val(A) -> .
{find the valuation of A and a reduced elt B with A~B.  Return B as an elt}
  if A`frame`is_first then
    return A;
  else
//"elt_val";
//"A`frame`segment",A`frame`segment;    
    phi_val := A`frame`segment`val;
//"phi_val",phi_val;
    coeff_with_vals := [elt_val(a) : a in A`coeffs];
//"coeff_with_vals",coeff_with_vals;
    this_val, pos := Minimum([phi_val*A`exponents[i]+coeff_with_vals[i]`valuation : 
                              i in [1..#A`exponents]]);
    return rec < elt_recformat() | 
                 frame := A`frame,
                 exponents := [A`exponents[pos]],
                 coeffs := [coeff_with_vals[pos]],
                 valuation := this_val 
               >;              
  end if;
end intrinsic;

intrinsic elt_to_poly(A,Ox) -> .
{}
// elt must be linear
//"elt_to_poly val",elt_val(A)`valuation;
  if A`frame`is_first then
    Pi := UniformizingElement(CoefficientRing(Ox));
    return A`unit*Pi^A`valuation;
  else
//Parent(A`frame`segment`frame`phi^A`exponents[1],"*";
//Parent(elt_to_poly(A`coeffs[1]));
//"Ox",Ox;
    return Ox!A`frame`segment`frame`phi^A`exponents[1]*Ox!elt_to_poly(A`coeffs[1],Ox);
  end if;
end intrinsic;

function constant_coefficient(polynomial)
  if polynomial eq 0 then
    return Parent(polynomial)!0;
  else
    return ConstantCoefficient(polynomial);
  end if;
end function;

intrinsic elt_from_poly(frame,f) -> .
{}
// return the (tree) elt representation of f with deg(f) < deg(frame`Phi)
  if frame`is_first then
    if f eq 0 then
      return rec< elt_recformat() | frame := frame, unit := 0, valuation := Infinity() >;
    else
      OK := CoefficientRing(Parent(frame`Phi));
      c := constant_coefficient(f);
//"Valuation(c)", Valuation(c);
      v := Integers()!Valuation(c);
      Pi := UniformizingElement(OK);
      return rec< elt_recformat() | frame := frame, unit := c div Pi^v, valuation := v >;
    end if;
  else
    expansion := phi_expansion(f,frame`segment`frame`phi);
    return rec< elt_recformat() | 
                frame := frame, 
                exponents := [0..#expansion-1],
                coeffs := [elt_from_poly(prev_frame(frame),b) : b in expansion]
              >;
  end if;
end intrinsic;

intrinsic elt_pow_raw(A,n) -> .
{}
// only for linear elements, that is coeffs and exps list have length 1
  if A`frame`is_first then
    return rec< elt_recformat() | 
                frame := A`frame, 
                unit := A`unit^Integers()!n,
                valuation := A`valuation*n  
              >;
  else
    return rec < elt_recformat() | 
                frame := A`frame, 
                exponents := [A`exponents[1]*n], 
                coeffs := [elt_pow_raw(A`coeffs[1],n)]
                //valuation := A`valuation*n  
              >;
   end if;

end intrinsic;

intrinsic elt_mult_raw(A,B) -> .
{}
// only for linear elements, that is coeffs and exps list have length 1
//"\nelt_mult_raw";
//"A";elt_print(A);
//"\nB";elt_print(B);
  if A`frame`depth ne B`frame`depth then
    error "OM elt_mult_raw: Elements must live in the same frame";
  end if;
  if A`frame`is_first then
    return rec< elt_recformat() | 
                frame := A`frame, 
                unit := A`unit*B`unit,
                valuation := A`valuation+B`valuation
                >;
  else
    return rec < elt_recformat() | 
                frame := A`frame, 
                exponents := [A`exponents[1]+B`exponents[1]], 
                coeffs := [elt_mult_raw(A`coeffs[1],B`coeffs[1])]
                //valuation := A`valuation+B`valuation
              >;
   end if;
end intrinsic;

intrinsic elt_mult(A,B) -> .
{A*B}
// only for linear elements, that is coeffs and exps list have length 1
// reduces exponents  
//"\nelt_mult";
 if A`frame`is_first then
    return elt_mult_raw(A,B);
  else
    exp := A`exponents[1]+B`exponents[1];
    if exp ge Degree(A`frame`phi) or exp lt 0 then
      // the exponent is not reduced
//"reduce !";
//"A";elt_print(A);
//"\nB";elt_print(B);
      s, exp := Quotrem(Integers()!exp, A`frame`segment`Eplus);
//"\nEplus",A`frame`segment`Eplus;
      coeff := elt_mult_raw(A`coeffs[1],B`coeffs[1]);
//"\ncoeff";elt_print(coeff);
//"\ngamma_psi";elt_print(A`frame`gamma_psi);
//"s",s;
      coeff := elt_mult(coeff,elt_pow_raw(A`frame`gamma_psi,s));
//"\ncoeff";elt_print(coeff);
    else  
      coeff := elt_mult(A`coeffs[1],B`coeffs[1]);
//"\ncoeff";elt_print(coeff);
    end if;
    return rec < elt_recformat() | 
                frame := A`frame, 
                exponents := [exp], 
                coeffs := [coeff] 
              >;
  end if;
end intrinsic;

intrinsic elt_residue(A) -> .
{A mod p as an element of O}
// only for linear OM elts
//"elt_residue(A)",A;
  if A`frame`is_first then 
    if A`valuation ne 0 then
      return 0;
    else
      return A`unit;
    end if;
  else
    if A`exponents[1] gt 0 then
      return 0;
    else
      return elt_residue(A`coeffs[1]);
    end if;
  end if;
end intrinsic;

intrinsic elt_mult_unit(A,gamma) -> .
{}
// multiply A by the scalar unit gamma
  if A`frame`is_first then
    B := A;
    B`unit := A`unit*gamma;
  else
    B := A;
    B`coeffs := [elt_mult_unit(a,gamma) : a in A`coeffs ];
  end if;
  return B;
end intrinsic;

////////////////////////////////////
// frame
///////////////////////////////////

intrinsic frame_recformat() -> .
{} 
   return recformat
        < is_first,        // true if this is the first frame
          depth,           // depth of the frame, first frame has depth 1
          segment,         // reference back to segment
          OK,              // base ring
          Phi,             // polynomial to be factored
          phi,             // approximation of factor that is under investigation
          residual_factor, // factor of residual polynomial
          gamma,           // lift of root of factor
          gamma_psi,       // gamma*Psi as an OM elt
          complete,        // true if this frame corresponds to a 
                           // unique irreducible factor of Phi
          F,               // maximum known inertia degree
          Fplus,           // increase in maximum known inertia degree
          degree           // maximum degree of irreducible factor 
                           // associated to the residual_factor
         >;
end intrinsic;

intrinsic OmFirstFrame(Phi) -> .
{}
  Kx<x> := Parent(Phi);

  return rec < frame_recformat() |
               is_first := true,
               phi := x,
               depth := 1,
               Phi := Phi,
               F := 1
               >;
end intrinsic;



intrinsic prev_frame(frame) -> .
{}
  return frame`segment`frame; 
end intrinsic;



////////////////////////////////////
// segment
////////////////////////////////////

intrinsic segment_recformat() -> .
{}
  return recformat
         < frame,    // reference back to frame
           slope,    
           val,      // -slope
           E,        // maximum known ramification
           Eplus, 
           degree,   // maximum degree of irreducible 
                     // factor associated with this segment
                     // length = (deg(frame`phi))*(length of segment)
           Psi,      // element with valuation Eplus*val
           //phi_val,
           residual_polynomial
           >;
end intrinsic;

////////////////////////////////////

intrinsic OmSegments(frame) -> .
{Return the list of segments of the newton polygon of phi}
vprint OM,1:"------------- OmSegments"; 
  Phi := frame`Phi;
  phi := frame`phi;
vprint OM,2:"phi",phi;
  OKy<y> := Parent(Phi);
  OK := CoefficientRing(OKy);
  RK, m := ResidueClassField(OK);
  RKz<z> := PolynomialRing(RK);
  mm := hom<OKy -> RKz | m, z>;

  if frame`is_first then
    E := 1;
  else
vprint OM,3:"frame",frame;
    E := frame`segment`E;
  end if;

  if frame`is_first then
    Phi_expanded := phi_expansion(Phi,phi);
  else  
//"frame`degree",frame`degree;
//"deg(phi)",Degree(phi);
//"phi",phi;
    //Phi_expanded := phi_expansion(Phi,phi:length:=frame`segment`degree/(Degree(phi)*frame`F)+1);
    Phi_expanded := phi_expansion(Phi,phi:length:=frame`degree/Degree(phi)+1);
  end if;
//"Phi_expanded",Phi_expanded;
  coeff_elts := [elt_from_poly(frame,g): g in Phi_expanded];
//"coeff_elts",coeff_elts;  
  A := [elt_val(a): a in coeff_elts];
  points := [<i,Rationals()!A[i]`valuation>: i in [1..#A] | A[i]`valuation ne Infinity()];
  // avoid points with infinite y-coordinate
//"points",points;Parent(points);  
  newton_polygon := NewtonPolygon(points);
  vertices := LowerVertices(newton_polygon);
  slopes := Slopes(newton_polygon)[1..#vertices-1];

//"vertices",vertices, ExtendedType(vertices);  
vprint OM,3:"segments: slopes",[s: s in slopes];

  if points[1][1] ne 1 then 
    // the first segment has slope -infinity
vprint OM,4:"slope",-Infinity();
    seg := rec < segment_recformat() | frame := frame >;
    if points[1][1] ne 2 then
      error "OmSegments: Segment with slope -infinity and length",points[1][1]+1, "which is greater than one.";
    end if;
    seg`degree := Degree(phi); 
    seg`slope := -Infinity();
    seg`Eplus := 1;
    seg`E := E;
    seg`val := Infinity();
    seg`residual_polynomial := z;
    segments := [seg];
  else
    // no segment with slope -infinity
    segments := [];
  end if;
  for i in [1..#slopes] do
    seg := rec < segment_recformat() | frame := frame >;
    seg`slope := slopes[i];
//"slope",seg`slope;
    if not frame`is_first then
      if frame`segment`slope le seg`slope then
        error "OmSegments -- The approximation has not improved.";
      end if;
    end if;
    seg`val := -seg`slope;
    e := Denominator(seg`slope);
    seg`Eplus := Integers()!(e/GCD(E,e)); 
    seg`E := E*seg`Eplus;
//"F",frame`F;
//"Eplus",seg`Eplus;
//"E",seg`E;
    seg`Psi := elt_with_val(frame,seg`Eplus*seg`val);
//"val",seg`Eplus*seg`val;
//"Psi val",elt_val(seg`Psi)`valuation;
//print "Psi";elt_print(seg`Psi);
    left_x := Integers()!vertices[i][1];
    right_x := Integers()!vertices[i+1][1];
//"right_x",right_x;
//Parent(right_x);
//"A",A;
//A[right_x];
    len_x := right_x-left_x;
//"len_x",len_x;
//"deg phi",Degree(phi);
    seg`degree := len_x*Degree(phi);
//"degree",seg`degree;
//seg`degree := 0;
    len := len_x/seg`Eplus;
    A_right_inv := elt_pow_raw(A[right_x],-1);
    R := OKy![ elt_residue(
                 elt_mult(A_right_inv,
                   elt_mult_raw(A[left_x+i*seg`Eplus],
                     elt_pow_raw(seg`Psi,i-len)))):
                       i in [0..len]]; 
//"R",R;
//"Parent(R)",Parent(R);
//"RKz",RKz;
    //seg`residual_polynomial := RKz!R;
//"mm",mm;    
    seg`residual_polynomial := mm(R);
    Append(~segments,seg);
  end for;
vprint OM,2:"--------- End Segments"; 
  return segments;

end intrinsic;
  
//////////////////////////////////////////////////////
// residual_factor
/////////////////////////////////////////////////////


intrinsic OmNextFrames(segment) -> .
{the next frames}
vprint OM,1:"------------- Frames"; 
  OK := CoefficientRing(Parent(segment`frame`Phi));
  OKx<x> := PolynomialRing(OK);
  RK, rm := ResidueClassField(OK);
  
vprint OM,2:"next frames: slope",segment`slope;  
  if segment`slope eq -Infinity() then
  // segment has slope -infinity, so phi is a factor of Phi
    frame := rec < frame_recformat() | is_first := false >;
    frame`degree := segment`degree;
    frame`phi := segment`frame`phi;
    frame`Phi := segment`frame`Phi;
    frame`Fplus := 1;
    frame`F := segment`frame`F;
    frame`segment := segment;
vprint OM,3:"--------- End Frames"; 
    return [frame];
  end if;
Parent(segment`residual_polynomial);
  factors := Factorization(segment`residual_polynomial);
  
  frames := [];
vprint OM,3:"old phi",segment`frame`phi; 
//"seg degree",segment`degree;
vprint OM,3:"rho factors",factors;
  for rho in factors do
//"  --- rho",rho;
    frame := rec < frame_recformat() | is_first := false, residual_factor := rho[1] >;
    frame`segment := segment; // this might change when new phi is computed
    frame`Fplus := Degree(rho[1]);
vprint OM,3:"deg rho",Degree(rho[1]);
//"mul rho",rho[2];
    frame`F := segment`frame`F * frame`Fplus;
    frame`degree := rho[2]*Degree(segment`frame`phi)*segment`Eplus;
vprint OM,3:"  --- frame degree",frame`degree;
    if frame`degree gt segment`degree then 
      error "OmFrames -- degree inconsistency"; 
    end if;   
    Append(~frames,frame);
  end for;
vprint OM,3:"--------- End Frames"; 
  return frames;
end intrinsic;


intrinsic om_next_phi(frame) -> .
{}  
vprint OM,1:"--------- Next phi";
  segment := frame`segment; 
  OK := CoefficientRing(Parent(segment`frame`Phi));
  OKx<x> := PolynomialRing(OK);
  RK, rm := ResidueClassField(OK);
//"frame`segment`slope",frame`segment`slope;
  rho := frame`residual_factor;
  if frame`Fplus eq 1 then
      frame`gamma := ConstantCoefficient(rho)@@rm;
vprint OM,3:"Gamma", frame`gamma; 
      frame`Phi := segment`frame`Phi;
      gamma_psi := elt_mult_unit(segment`Psi,frame`gamma);
//"val Psi",elt_val(segment`Psi)`valuation;      
//vprint OM,3:"\nPsi";elt_print(segment`Psi); 
//vprint OM,4:"\ngamma_psi";elt_print(gamma_psi); 
//vprint OM,4:"\ngamma_psi\n",elt_to_poly(gamma_psi,Parent(segment`frame`phi));
      frame`phi := segment`frame`phi^segment`Eplus+elt_to_poly(gamma_psi,Parent(segment`frame`phi));
vprint OM,2:"  --- new phi",frame`phi; 
    else
      OL := CyclotomicUnramifiedExtension(OK,Degree(rho));
      OLY<Y> := PolynomialRing(OL);
//"OL",OL;      
      frame`gamma := OL.1;
//"\nGamma";frame`gamma; 
      frame`OK := OL;
      gamma_psi := elt_mult_unit(segment`Psi,frame`gamma);
//"\ngamma_psi";elt_print(gamma_psi); 
      frame`Phi := OLY!segment`frame`Phi;
      frame`phi := OLY!segment`frame`phi^segment`Eplus-elt_to_poly(gamma_psi,OLY);
vprint OM,2:"  --- new phi",frame`phi; 
    end if;
    if segment`Eplus eq 1 and not segment`frame`is_first then
      // this was an improvement step, skip the previous frame
vprint OM,2:"improvement step";
      frame`segment := segment`frame`segment;
      frame`depth := segment`frame`depth;
      // gamma_psi depends on segment of previous frame, save it
      frame`gamma_psi := segment`frame`gamma_psi;
//    elif segment`Eplus eq 1 and segment`frame`is_first then
      // this was an improvement step, the previous frame was the first frame
      // the is the first frame now.
//      frame`is_first := true;
    else
      //frame`segment := segment;
      frame`depth := segment`frame`depth+1;
      frame`gamma_psi := gamma_psi;
    end if;
vprint OM,4:"depth",frame`depth;
    
    return frame;

end intrinsic;


//////////////////////////////////////////////////////
// main
/////////////////////////////////////////////////////


intrinsic om_tree_sub(frame) -> .
{}
  segs := OmSegments(frame);
//"#segs",#segs;
  if segs[1]`slope eq -Infinity() then
    F := OmNextFrames(segs[1]);
  else
    F := [];
  end if;
  new_frames := [OmNextFrames(s) : s in segs | s`slope ne -Infinity()];
//"#frames",[#a: a in new_frames];
  new_frames := [om_next_phi(f): f in Flat(new_frames)];
//"om tree sub #frames",#new_frames;
  for i in [1..#new_frames] do
   f := new_frames[i];
//" FFFFFFFFFFFFFFFFFFFFFFFFFFF frame",i;//f;
//" slope",f`segment`slope;
//" E",f`segment`E;
//" Phi",f`Phi,Parent(f`Phi);
//" irred factor degree bound",f`degree;
//" deg phi",Degree(f`phi);
//" phi",f`phi;
    if f`degree eq Degree(f`phi) then
//" DONE !";
      // f corresponds to an irreducible factor
      Append(~F,f);
    else
//" CONTINUE";
      F cat:= om_tree_sub(f);
//"#F",#F;
    end if;
  end for;
  return F;
end intrinsic;


intrinsic OmTree(Phi) -> .
{}
  first_frame := OmFirstFrame(Phi);
  return om_tree_sub(first_frame);
end intrinsic;

///////////////////////////////////////////////////
// Lifting and Single Factor
///////////////////////////////////////////////////

function om_index_sub(frame)
  if frame`is_first then
    return 1;
  else
    return frame`segment`Eplus*frame`segment`slope+om_index_sub(frame`segment`frame);
  end if;
end function;

intrinsic OmIndex(frame) -> .
{}
  return Floor(om_index_sub(frame));
end intrinsic;

intrinsic OmLift(frame,nu) -> .
{Lift an approximation to an irreducible factor, to a factor with precision nu}
  // ??? fix precision
//"lifting";
//"frame`phi",frame`phi;
 
  if frame`segment`val eq Infinity() then
  // the segment of phi has slope -infinity, so phi is a factor already.
    return frame`phi; 
  end if;
  
  prec := nu+4*OmIndex(frame)+Ceiling(frame`segment`slope);; 

  Ox<x> := Parent(frame`phi);
  O     := CoefficientRing(Ox);
  piO   := UniformizingElement(O);
  
//U := BaseRing(O); 
 
  Q := O;
  Qz := Ox;
  //Q     := quo<O | piO^prec>;
  //Qz<z> := PolynomialRing(Q);

  pi := UniformizingElement(Q);
//"pi",pi;
  phi := Qz!frame`phi;
//"phi",phi;
  Phi := Qz!frame`Phi;

  phi := frame`phi;
  a, a0 := Quotrem(frame`Phi,phi);
//elt_a0 := elt_val(elt_from_poly(frame,a0));
//"val(a0)",elt_a0`valuation; 
  a1 := a mod phi;
//"a1  ",a1;
//"a1 elt"; elt_print(elt_from_poly(frame,a1));"\n";
//"a1<>",elt_to_poly(elt_from_poly(frame,a1),Ox);
  elt_a1 := elt_val(elt_from_poly(frame,a1));
//"val(a1)",elt_a1`valuation; 
//"a1 as elt";elt_print(elt_a1);"\n";
  elt_psi, exp := elt_with_val(frame,-elt_a1`valuation:extval:=true);
//"elt_psi";elt_print(elt_psi);"\n";
//"exp",exp;
  psi := elt_to_poly(elt_psi,Ox);
//"psi",psi;
  A0 := ((psi*a0) mod phi) div (pi^-exp);
//"a1<v>",elt_to_poly(elt_a1,Ox);
  A1 := ((psi*elt_to_poly(elt_a1,Ox)) mod phi) div (pi^-exp);
//"A1",A1;
//"psi*a1 mod phi div pi^-exp",Qz![a div pi^-exp:a in Eltseq((psi*a1) mod phi)];;
//"psi*a1 mod phi",((psi*a1) mod phi);;
//"psi*a1 mod phi",(psi*a1) mod phi;
//"psi*a1 mod phi div pi^-exp",Qz![a div pi^-exp:a in Eltseq((psi*a1) mod phi)];;
//"v(a1)", elt_val(elt_a1)`valuation;
  A1_elt := elt_from_poly(frame,A1);
//"v(A1)", elt_val(A1_elt)`valuation;
//"C1",elt_residue(A1_elt);
  C1_inv := O!elt_residue(A1_elt)^-1;
//"A0",A0;
//"C1_inv",C1_inv;
  phi := phi +((A0*C1_inv) mod phi);
  h := 1;
  while h le frame`segment`E*nu do
    c, c0 := Quotrem(frame`Phi,phi);
    c1 := c mod phi;
    C0 := (psi*c0 mod phi) div pi^-exp;
    C1 := (psi*c1 mod phi) div pi^-exp;
    C1_inv := (C1_inv*(2-C1*C1_inv)) mod phi; 
    phi := phi + ((C0*C1_inv) mod phi);
    h := 2*h;
  end while;
//"phi",phi; 
  return phi;
end intrinsic;

intrinsic OmSingleFactor(Phi) -> .
{}
//"single factor";
  frame := OmFirstFrame(Phi);
  segs := OmSegments(frame);
  frame := OmNextFrames(segs[1])[1];
  while frame`degree ne Degree(frame`phi) do
    segs := OmSegments(frame);
    frame := OmNextFrames(segs[1])[1];
  end while;
  return OmLift(frame,Precision(CoefficientRing(Parent(Phi))));
end intrinsic;


//////////////////////////////////////////////////////////////////
// Characteristic Polynomial
//////////////////////////////////////////////////////////////////
// from RoundFour
//////////////////////////////////////////////////////////////////
// compute the characteristic polynomials using newton relations
// (cohen 1, p 161)
// the traces_from_poly of tyhe reference polynomials Phi are 
// precomputed and applied in th computation of char polys wrt
// Phi.        


traces_from_poly := function(Phi)

  n := Degree(Phi);
  R := CoefficientRing(Parent(Phi));

  // precision loss occurs in poly_from_traces by division
  if n gt 0 then
    inc := Maximum([Valuation(R!a):a in [1..n]]);
    inc := inc^2;
  else
    inc := 0;
  end if;
  //print "inc",inc;

  tmp_R := ChangePrecision(R,Precision(R)+inc);
  tmp_PR := PolynomialRing(tmp_R);
  tmp_Phi := tmp_PR!Phi;

  t := function(i)
    if i lt 0 then
      return 0;
    else
      return Coefficient(tmp_Phi,i);
    end if;
  end function;

  S := [];

  for k := 1 to n do
    Sk := -k*t(n-k);
    for i := 1 to Minimum(n,k-1) do
      Sk := Sk - t(n-i)*S[k-i];
    end for;
    Append(~S,Sk);
  end for;
  return S;
end function;

poly_from_traces := function(S)

  n := #S;
  Tn := [];

  for k := 1 to n do
    Tnk := -S[k];
    for i := 1 to k-1 do
      Tnk := Tnk - Tn[i]*S[k-i];
    end for;
    //print Tnk;
    //print "mod", Tnk mod k;
    Tnk := Tnk div k;
    //print k,Tnk;
    Append(~Tn,Tnk);
  end for;
  Tn := Reverse(Tn);
  Append(~Tn,Parent(S[1])!1);
  return PolynomialRing(Parent(S[1]))!Tn;

end function;

oldchi := function(Phi,theta)
// The char poly of theta wrt Phi -- Kernel version
  PR := Parent(Phi);
  R := CoefficientRing(PR);
  //print "chi Phi", Phi;
  //print "chi theta",theta;

  if theta eq PR.1 then
    return Phi;
  end if;

  ML := [];
  seq := [ R!0 : a in [1..Degree(Phi)]];
  seq[1] := R!1;
  Append(~ML, seq);
  poly := theta mod Phi;


  for i in [1..Degree(Phi)] do
    //print i;
    if Degree(poly) lt Degree(Phi)-1 then
      seq := Eltseq(poly+PR.1^(Degree(Phi)-1));
      seq[Degree(Phi)] := seq[Degree(Phi)]-1;
    else
      seq := Eltseq(poly);
    end if;
    seq := seq[1..Degree(Phi)];
    Append(~ML, seq);
    M := Matrix(Reverse(ML));
    //print "M",M;
    K := Kernel(M);
    //print "K",K;
    if #Basis(K) gt 0 then
      //print "K",K;
      res := PR![e : e in Reverse(Eltseq(Basis(K)[1]))];
      //print res;
      //print Degree(res);
      //print LeadingCoefficient(res);
      //print res;
      //if Degree(Phi) mod Degree(res) eq 0 then
         //print "res1",res;
        if IsUnit(LeadingCoefficient(res)) then
           res := res div LeadingCoefficient(res);
           //print "res2",res;
           //print "basislen",#Basis(K);
           //print K;
           //print "chi kernel",res;
           return res;
        end if;
      //end if;
    end if;
    poly := poly*theta mod Phi;
  end for;
  error "chi";
end function;

chi := function(Phi, traces_Phi, theta)

  if Characteristic(CoefficientRing(Phi)) ne 0 then
    return oldchi(Phi,theta);
  end if;

  n := Degree(Phi);
  assert n gt 0;

  tmp_PR := PolynomialRing(Universe(traces_Phi));
//"tmp_PR",tmp_PR;
//"theta",Parent(theta);  
  tmp_theta := tmp_PR!theta;
  tmp_Phi := tmp_PR!Phi;

  S := [];

  pow := 1;
  for i := 1 to n do
    pow := pow*tmp_theta mod tmp_Phi;
    //print "pow",pow;
    Si := 0;
    Si := n*Coefficient(pow,0);
    for j := 1 to Degree(pow) do
      if Coefficient(pow,j) ne 0 then
        Si +:= traces_Phi[j]*Coefficient(pow,j);
      end if;
    end for;
    Append(~S,Si);
  end for;
  //print S;
  //return S;
  return Parent(Phi)!poly_from_traces(S);
end function;


chi_with_den := function(Phi, traces_Phi, theta, denval)
//{Find char poly of theta * pi^val}
//"chi";
    PR := Parent(Phi);
    R := CoefficientRing(Phi);
//"R",R;
    tmp_R := ChangePrecision(R,Precision(R)+4*denval*Degree(Phi));
    tmp_PR := PolynomialRing(tmp_R);
    res := chi(Phi, traces_Phi, theta);
    res := tmp_PR!res;
    pi := UniformizingElement(tmp_R);
    res := Evaluate(res, pi^denval * tmp_PR.1);
    res div:= LeadingCoefficient(res);
    return res;
    assert IsMonic(res);
    return PR!res;
end function;



//////////////////////////////////////////////////////////
// Extensions
/////////////////////////////////////////////////////////


intrinsic OmDefiningPolynomial(frame) -> .
{}
  Ox<x> := Parent(frame`phi);
  O := CoefficientRing(Ox); 
  
  prec := Maximum(Precision(O),1+2*Valuation(O!Degree(frame`phi)));
  f := OmLift(frame,prec);
  
  elt_psi, psi_pi_exp := elt_with_val(frame,1/frame`segment`E : extval:=true);
  psi := Parent(f)!elt_to_poly(elt_psi,Ox);
//"psi",psi,"* pi ^",psi_pi_exp; 
 
  traces_f := traces_from_poly(f);
  this_chi := chi_with_den(f, traces_f, psi, -psi_pi_exp); 
  return this_chi;

end intrinsic;
/////////////////////////////////////////////////////////////////
// We assume that the basefield is a wildliy ramified extension
// over a tamely ramified extension
// over an unramified extension
/////////////////////////////////////////////////////////////////

function change_elt_precision(a,prec)
  o := Parent(a);
  oo := ChangePrecision(o,prec);
  return o!oo!a; 
end function;

intrinsic reduce_poly(f) -> .
{}
// f must be eisenstein
//"reduce_poly in",f;
  o := CoefficientRing(Parent(f));
  n := Degree(f);
  fseq := Eltseq(f);
  j := Min([n*Valuation(o!n)] cat [n*(Valuation(o!i)+Valuation(fseq[i+1])-1)+i:i in [1..n-1]]);
  prec := Ceiling((n+2*j+1)/n);
  oo := ChangePrecision(o,prec);
//"reduce_poly out",Polynomial(o,Polynomial(oo,f)); 
  return Polynomial(o,Polynomial(oo,f)); 
end intrinsic;

function map1(x,a)
  l := Eltseq(x);
  return &+[a^(i-1)*l[i]:i in [1..#l]];
end function; 

function map2(x,a,b)
  l := Eltseq(x);
  return &+[b^(i-1)*map1(l[i],a):i in [1..#l]];
end function; 




intrinsic OptimizedUnramifiedExtension(R,F) -> .
{unramified extension of degree F of R}
//"optimized unramified extension";
  if Characteristic(R) ne 0 then
    P := CoefficientRing(R);
    K := ext<P|F>;
    return PowerSeriesRing(K);
  end if;


  P := PrimeRing(R);
  f := InertiaDegree(R,P);
  // e := RamificationIndex(R,P);
  
  if f eq 1 then 
    U := CyclotomicUnramifiedExtension(P,F);
    if P cmpeq R then
      return U;
    elif BaseRing(R) cmpeq P then
      // R is a ramified extension
      phi := DefiningPolynomial(R);
//"phi",phi;      
      return TotallyRamifiedExtension(U,reduce_poly(Polynomial(U,phi)));
    else
      // R is a wildly ramified extension over a tamely ramified extension
      R_tame := BaseRing(R);
      phi := DefiningPolynomial(R);
      psi := DefiningPolynomial(R_tame);
      T := TotallyRamifiedExtension(U,reduce_poly(Polynomial(U,psi)));
      new_phi := Polynomial(T,[T!Eltseq(a):a in Eltseq(phi)]);
      return TotallyRamifiedExtension(T,reduce_poly(new_phi));
    end if;
  else // f gt 1
    return CyclotomicUnramifiedExtension(P,f*F);
    if BaseRing(R) eq P then
      // R is unramified
      return CyclotomicUnramifiedExtension(P,f*F);
    elif BaseRing(BaseRing(R)) eq P then
      // R is ramified over unramified 
      R_un := BaseRing(R);
      U := CyclotomicUnramifiedExtension(P,f*F);
      rho := DefiningPolynomial(R_un);      
      a := Roots(rho:Max:=1)[1][1];
      phi := DefiningPolynomial(R);
      new_phi := Polynomial(U,[map1(b,a):b in Eltseq(phi)]);
      return TotallyRamifiedExtension(U,reduce_poly(new_phi));
    else // R is wildly ramified over tamely ramified over unramified 
      R_tame := BaseRing(R);
      R_un := BaseRing(R_tame);
      U := CyclotomicUnramifiedExtension(P,f*F);
      rho := DefiningPolynomial(R_un);      
      a := Roots(rho:Max:=1)[1][1];
      phi := DefiningPolynomial(R_tame);
      new_phi := Polynomial(U,[map1(c,a):c in Eltseq(phi)]);
      T := TotallyRamifiedExtension(U,reduce_poly(new_phi));
      psi := DefiningPolynomial(R);
      new_psi := Polynomial(T,[map2(c,a,T.1):c in Eltseq(phi)]);
      return TotallyRamifiedExtension(U,reduce_poly(new_psi));
    end if;
end if;

end intrinsic;


intrinsic OptimizedTamelyRamifiedExtension(R,E0,C) -> .
{tamely ramified extension of R given by x^E0+C}
//"optimized tamely ramified extension";
  
  P := PrimeRing(R);
  p := Prime(R);
  e := RamificationIndex(R,P);
  w := Valuation(e,p); 
  e0 := e div p^w;  
  pi := UniformizingElement(R); 

  if RamificationIndex(R,P) eq 1 then
    Rz<z> := PolynomialRing(R); 
    C:=Rz!C;
    return TotallyRamifiedExtension(R,reduce_poly(z^E0+C));
  end if;

  e := RamificationIndex(R,P);
  w := Valuation(e,p); 
  e0 := e div p^w;  

   if w eq 0 then
//"w=0";
    U := BaseRing(R);
    UC := Norm(C);
    Uz<z> := PolynomialRing(U);
    return TotallyRamifiedExtension(U,reduce_poly(z^(e0*E0)+UC));
  elif e0 eq 1 then // e = p^w
//"e0=1";   
    U := BaseRing(R);
    UC := Norm(C);
    Uz<z> := PolynomialRing(U); 
    T := TotallyRamifiedExtension(U,reduce_poly(z^E0+U!UC));
    Ty<y> := PolynomialRing(T);
    phi := Ty!DefiningPolynomial(R);
    one, a, b := XGCD(E0,p^w);
    if a lt 1 then
      a := a+p^w; b := b-E0;
    end if;
    PI := y^a; piden := -b;
    traces_phi := traces_from_poly(phi);
    chi := chi_with_den(phi, traces_phi, PI, piden); 
    return TotallyRamifiedExtension(T,reduce_poly(chi));
  else // e = e0 * p^w
//"else";   
    
    U := BaseRing(BaseRing(R));
//"U",U;
    Uy<y> := PolynomialRing(U); 

    pi := UniformizingElement(R);
//"R!C div pi",(R!C div pi);
    epsilon := U!Eltseq(Eltseq(R!C div pi)[1])[1];
//"epsilon",Parent(epsilon); 

    R_tame := BaseRing(R);
    phi := DefiningPolynomial(R); // in BaseRing(R);
    pi0 := UniformizingElement(R_tame);
    gamma0 := U!Eltseq(ConstantCoefficient(phi) div pi0)[1];
//"gamma0",Parent(gamma0); 
    
    phi0 := DefiningPolynomial(R_tame);
    delta := ConstantCoefficient(phi0);
//"delta",Parent(delta); 
    Delta := delta*gamma0^e0*epsilon^e;
//"Delta",Parent(Delta);
    T := TotallyRamifiedExtension(U,reduce_poly(y^(E0*e0)+Delta));
    Tz<z> := PolynomialRing(T); 
//"phi",phi;
    Tpi0 := T.1^E0 div (epsilon^(p^w)*gamma0);
    new_phi := Polynomial(T,[map1(c,Tpi0):c in Eltseq(phi)]);
//"new_phi",new_phi;

    one, a, b := XGCD(E0,p^w);
    if a lt 1 then
      a := a+p^w; b := b-E0;
    end if;
//"E0",E0;
//"p^w",p^w;
//"one",one;
//"a",a;
//"b",b;   
    PI := z^a; piden := -b;
    traces_new_phi := traces_from_poly(new_phi);
    chi := chi_with_den(new_phi, traces_new_phi, PI, piden); 
    return TotallyRamifiedExtension(T,reduce_poly(chi));
  end if;

end intrinsic;



intrinsic OptimizedWildlyRamifiedExtension(R,PHI) -> .
{An optimized representation of the wildly ramified extension of R generated by PHI, i.e., as a wildly ramified extension over a tamely ramified extension over an unramifed extension of Zp.}
//"optimized wildly ramified extension";
Ry<y>:=Parent(PHI);
//"PHI",PHI;  
  P := PrimeRing(R);
  p := Prime(R);
  e := RamificationIndex(R,P);
  w := Valuation(e,p); 
  e0 := e div p^w;  
 
  if e eq 1 or  w eq 0 then
    // R is unramified or R is tamely ramified
vprint OM,3:"  no previous wilderness";
    return TotallyRamifiedExtension(R,reduce_poly(PHI));
  else 
    // R is wildly ramified
    T := BaseRing(R);
vprint OM,3:"  reducing";
   PHI := reduce_poly(PHI);
//"  PHI",PHI,"in",Parent(PHI)," Eisenstein ?",IsEisenstein(PHI);
vprint OM,3:"  norm start";
//    PSI := Norm(PHI);
   PSI := DefiningPolynomial(TotallyRamifiedExtension(R,PHI),T);
vprint OM,3:"  norm done";    
    // PSI is in T[y]
    PSI:= reduce_poly(PSI);
//"  PSI",PSI;
    return TotallyRamifiedExtension(T,PSI);
  end if;

end intrinsic;


////////////////////////////////////////////////////////////////////////
// Splitting Field
////////////////////////////////////////////////////////////////////////


intrinsic om_splitting_sub(frame,Phi,bound) -> .
{}
//"Parent Phi",Parent(Phi);
// "split: segs";
  segs := OmSegments(frame);
// "split: frames";
  new_frames := [OmNextFrames(s) : s in segs];
  new_frames := Flat(new_frames);
//  new_frames := PermuteSequence(new_frames,Random(SymmetricGroup(#new_frames)));
 
  // new inertia ?
  F := LCM([f`Fplus:f in new_frames]);
  if F gt 1 then
    K := OptimizedUnramifiedExtension(CoefficientRing(Parent(frame`phi)),F);
    return om_splitting_start(Polynomial(K,Phi),Phi,bound);
  end if;
      
  for new_frame in new_frames do
k := CoefficientRing(Parent(new_frame`segment`frame`Phi)); pk := PrimeRing(k); 
//"k",k;
vprint OM,3:"e_K :",RamificationIndex(k,pk)," f_K :",InertiaDegree(k,pk)," factor degree :",new_frame`degree," slope :",new_frame`segment`slope," E :",new_frame`segment`E," F :",new_frame`F;    
    if new_frame`degree eq 1 then
vprint OM,3:"split: discard";        
    else
      this_frame := om_next_phi(new_frame);
vprint OM,3:" --  deg(phi) :",Degree(this_frame`phi),"phi",this_frame`phi;;    
      if this_frame`degree eq Degree(this_frame`phi) then
//"if";
        // this_frame corresponds to an irreducible factor
        if this_frame`degree gt 1 then
//"ifif";
          // the irreducible factor is not linear
vprint OM,2:"split: increase E";
          g := OmDefiningPolynomial(this_frame); 
          //p := Prime(CoefficientRing(Parent(this_frame`phi)));
          R := CoefficientRing(Parent(this_frame`phi));
	   state:=RamificationIndex(k,pk)*InertiaDegree(k,pk)*Degree(g);
	   "Degree of splitting field is a least", state,".";
          if bound ne 0 then
            if RamificationIndex(k,pk)*InertiaDegree(k,pk)*Degree(g) gt bound then
              return TotallyRamifiedExtension(R,g);
            end if;
          end if;
          p := Prime(CoefficientRing(Parent(this_frame`phi)));
          R := CoefficientRing(Parent(this_frame`phi));
          w := Valuation(Degree(g),p);
          E0 := Degree(g) div p^w;
          "optimizing extension";
          if w eq 0 then
             vprint OM,3:"tame";
             K := OptimizedTamelyRamifiedExtension(R,E0,ConstantCoefficient(g));
          elif E0 gt 1 then
             vprint OM,3:"tame";
            K := OptimizedTamelyRamifiedExtension(R,E0,ConstantCoefficient(g));
          else
            vprint OM,3:"wild";
            K := OptimizedWildlyRamifiedExtension(R,g);
          end if;
          //K := OptimizedRamifiedExtension(g);
//"K",K;
//"Parent Phi",Parent(Phi);
          "starting over new field",K;
          return om_splitting_start(Polynomial(K,Phi),Phi,bound);
        end if;
      else // this_frame`degree ne Degree(this_frame`phi)
vprint OM,2:"split: continue";
        K := om_splitting_sub(this_frame,Phi,bound);
        if K cmpne CoefficientRing(Parent(this_frame`phi)) then
          // we have the splitting field, get out of here
vprint OM,2:"split: done";
          return K;
        end if;
      end if;
    end if;
  end for;
vprint OM,1:"split: fallthrough";
  K := CoefficientRing(Parent(frame`phi));
  return K;
end intrinsic;

intrinsic om_splitting_start(Phi_now, Phi_base,bound) -> .
{}
vprint OM,1:"split: start";
  first_frame := OmFirstFrame(Phi_now);
  return om_splitting_sub(first_frame,Phi_base,bound);
end intrinsic;

intrinsic OmSplittingField(Phi:bound:=0) -> .
{}
vprint OM,1:"OmSplittingField of",Phi;
  O := CoefficientRing(Parent(Phi));
  //"O",O;
  if Characteristic(O) eq 0 then
vprint OM,1:"characteristic 0";
    Q := quo < O | UniformizingElement(O)^Precision(O)>;
    Qt<t> := PolynomialRing(Q);
    return om_splitting_start(Qt!Phi,Qt!Phi,bound);
  else
vprint OM,1:"characteristic",Characteristic(O);
    return om_splitting_start(Phi,Phi,bound);
  end if;
end intrinsic;
