174,0
T,ModTupFld_FldPadExact,ModTupFldElt_FldPadExact,1,StrPadExact
A,ModTupFld_FldPadExact,7,base_field,degree,generic,outer_generators,generators,dimension,zero
T,ModTupFldElt_FldPadExact,-,1,PadExactElt
A,ModTupFldElt_FldPadExact,4,negation,component,eltseq,norm
A,FldPadExact,1,vector_space
S,VectorSpace,The full vector space over K of dimension n,0,2,0,0,0,0,0,0,0,148,,0,0,FldPadExact,,ModTupFld_FldPadExact,-38,-38,-38,-38,-38
S,KSpace,The full vector space over K of dimension n,0,2,0,0,0,0,0,0,0,148,,0,0,FldPadExact,,ModTupFld_FldPadExact,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,ModTupFld_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,ModTupFldElt_FldPadExact,,-38,-38,-38,-38,-38,-38
S,InterpolateEpochs,Interpolates between the given epochs,0,4,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,148,,0,0,ModTupFldElt_FldPadExact,,168,-38,-38,-38,-38,-38
S,Ngens,The number of generators of V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,NumberOfNames,The number of generators of V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Name,The ith generator of V,0,2,0,0,0,0,0,0,0,148,,0,0,ModTupFld_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,.,The ith generator of V,0,2,0,0,0,0,0,0,0,148,,0,0,ModTupFld_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,AssignNames,Assigns names to the generators of V,1,1,1,82,0,298,2,0,0,1,0,0,0,0,82,,1,1,ModTupFld_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Vector,The vector of length n over F defined by cs,0,3,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Vector,The vector whose coefficients are given by cs,1,0,1,82,0,FldPadExactElt,1,0,0,0,0,0,0,0,82,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Zero,The zero vector,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ZeroVector,The zero vector,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ZeroVector,The zero vector,0,2,0,0,0,0,0,0,0,148,,0,0,FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of V. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,-1,,0,0,ModTupFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of V. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,ModTupFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of V. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,82,,0,0,ModTupFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,BaseField,The base field of V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,BaseField,The base field of v,0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,Degree,"If V is a subspace of `K^n`, returns `n`. That is, the number of columns in vectors in V",0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Generic,The generic vector space containing V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,ModTupFld_FldPadExact,-38,-38,-38,-38,-38
S,IsGeneric,"True if V is generic, i.e. it was created as the full-dimensional vector space with default generators",0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,36,-38,-38,-38,-38,-38
S,Dimension,The dimension of V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Generators,The generators of V,0,1,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,82,-38,-38,-38,-38,-38
S,ExistsCoveringStructure,True if there is a space containing both V1 and V2,0,2,0,0,0,0,0,0,0,ModTupFld_FldPadExact,,0,0,ModTupFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,CanChangeRing,True if the base ring of v can be changed to F,0,2,0,0,0,0,0,0,0,FldPadExact,,0,0,ModTupFldElt_FldPadExact,,36,ModTupFldElt_FldPadExact,-38,-38,-38,-38
S,ChangeRing,A copy of v with its base field changed to F,0,2,0,0,0,0,0,0,0,FldPadExact,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Nrows,The number of rows in v (always 1),0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,148,-38,-38,-38,-38,-38
S,Ncols,The number of columns in v,0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,148,-38,-38,-38,-38,-38
S,Eltseq,The components of v as a sequence,0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,Component,The ith component of v,0,2,0,0,0,0,0,0,0,148,,0,0,ModTupFldElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,@,The ith component of v,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,148,,FldPadExactElt,-38,-38,-38,-38,-38
S,-,Negation,0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,-1,,0,0,ModTupFldElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,-1,,0,0,ModTupFldElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,&+,Sum,1,0,1,82,0,ModTupFldElt_FldPadExact,1,0,0,0,0,0,0,0,82,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Scalar multiplication,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Scalar multiplication,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,FldPadExactElt,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,/,Scalar division,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,InnerProduct,Inner product,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,InnerProduct,Inner product,0,2,0,0,0,0,0,0,0,-1,,0,0,ModTupFldElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,InnerProduct,Inner product,0,2,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,Norm,Norm,0,1,0,0,0,0,0,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
