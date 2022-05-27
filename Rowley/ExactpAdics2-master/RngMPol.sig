174,0
T,RngMPol_FldPadExact,RngMPolElt_FldPadExact,1,StrPadExact
A,RngMPol_FldPadExact,7,base_ring,rank,varnames,generator,generators,zero,one
T,RngMPolElt_FldPadExact,-,1,PadExactElt
A,RngMPolElt_FldPadExact,5,weak_monomials,weak_coefficients,is_monomial,exponents,negation
S,PolynomialRing,The polynomial ring of rank k over F,0,2,0,0,0,0,0,0,0,148,,0,0,FldPadExact,,RngMPol_FldPadExact,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,RngMPol_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,RngMPolElt_FldPadExact,,-38,-38,-38,-38,-38,-38
S,BaseRing,The base ring of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,BaseRing,The base ring of f,0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,Rank,The rank of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,148,-38,-38,-38,-38,-38
S,eq,Equality,0,2,0,0,0,0,0,0,0,RngMPol_FldPadExact,,0,0,RngMPol_FldPadExact,,36,-38,-38,-38,-38,-38
S,Generator,The ith generator of R,0,2,0,0,0,0,0,0,0,148,,0,0,RngMPol_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Generators,The generators of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,82,-38,-38,-38,-38,-38
S,Name,The ith generator of R,0,2,0,0,0,0,0,0,0,148,,0,0,RngMPol_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,.,The ith generator of R,0,2,0,0,0,0,0,0,0,148,,0,0,RngMPol_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,NumberOfNames,The number of generators of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,148,-38,-38,-38,-38,-38
S,AssignNames,Assigns names to the generators of R,1,1,1,82,0,298,2,0,0,1,0,0,0,0,82,,1,1,RngMPol_FldPadExact,,-38,-38,-38,-38,-38,-38
S,AssignNames,Assigns names to the generators of xR from R,0,2,0,0,1,0,0,0,0,RngMPol_FldPadExact,,1,1,64,,-38,-38,-38,-38,-38,-38
S,SetApproximationHook,Called by SetApproximation,0,3,0,0,1,0,0,1,1,64,,0,0,148,,0,0,RngMPol_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Zero,The zero of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,One,The one of R,0,1,0,0,0,0,0,0,0,RngMPol_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of R. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,-1,,0,0,RngMPol_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of R. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPol_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of R. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,63,,0,0,RngMPol_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of R. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,82,,0,0,RngMPol_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of R. If so, also returns the coerced element",1,1,1,82,0,303,2,0,0,0,0,0,0,0,82,,0,0,RngMPol_FldPadExact,,36,-1,-38,-38,-38,-38
S,MonomialCoefficient,The coefficient of monomial m in f,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPolElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,ExponentsCoefficient,The coefficient of exponent e in f,1,1,1,82,0,148,2,0,0,0,0,0,0,0,82,,0,0,RngMPolElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,Monomial,The monomial of R with exponents e,1,1,1,82,0,148,2,0,0,0,0,0,0,0,82,,0,0,RngMPol_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,IsDefinitelyMonomial,"True if f is definitely a monomial (i.e. has one term). If so, also returns its exponents",0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,36,82,-38,-38,-38,-38
S,Exponents,"The exponents of m, which must be a monomial",0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,-,Negate,0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Add,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Add,0,2,0,0,0,0,0,0,0,-1,,0,0,RngMPolElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,+,Add,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,-,Subtract,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,-,Subtract,0,2,0,0,0,0,0,0,0,-1,,0,0,RngMPolElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,-,Subtract,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,*,Multiply,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Multiply,0,2,0,0,0,0,0,0,0,-1,,0,0,RngMPolElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,*,Multiply,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,/,Divide by scalar,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,/,Divide by scalar,0,2,0,0,0,0,0,0,0,-1,,0,0,RngMPolElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,/,Divide by scalar,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,63,,-1,-38,-38,-38,-38,-38
S,&+,Sum,1,0,1,82,0,RngMPolElt_FldPadExact,1,0,0,0,0,0,0,0,82,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,&*,Product,1,0,1,82,0,RngMPolElt_FldPadExact,1,0,0,0,0,0,0,0,82,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Derivative,The mth derivative of f with respect to v,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Derivative,The mth derivative of f with respect to v,0,3,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,148,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Derivative,The derivative of f with respect to v,0,2,0,0,0,0,0,0,0,148,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Derivative,The derivative of f with respect to v,0,2,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,0,0,RngMPolElt_FldPadExact,,RngMPolElt_FldPadExact,-38,-38,-38,-38,-38
S,Evaluate,Evaluates `f(xs)`,1,1,1,82,0,FldPadExactElt,2,0,0,0,0,0,0,0,82,,0,0,RngMPolElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,Evaluate,Evaluates `f(xs)`,0,2,0,0,0,0,0,0,0,82,,0,0,RngMPolElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,Evaluate,Evaluates `f(xs)`,1,1,1,82,0,FldPadExactElt,2,0,0,0,0,0,0,0,82,,0,0,63,,-1,-38,-38,-38,-38,-38
S,WeakMonomials,The monomials of f. Some of the corresponding coefficients may be zero,0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,WeakCoefficients,The coefficients of f corresponding to `WeakMonomials(f)`,0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,WeakCoefficientsAndMonomials,The coefficients and monomials of f,0,1,0,0,0,0,0,0,0,RngMPolElt_FldPadExact,,82,82,-38,-38,-38,-38
S,IsHenselLiftable,"True if xs are Hensel liftable to roots of fs. If so, also returns the lifted roots",2,0,1,82,0,RngMPolElt_FldPadExact,1,1,82,0,FldPadExactElt,2,0,0,0,0,0,0,0,82,,0,0,82,,36,82,-38,-38,-38,-38
