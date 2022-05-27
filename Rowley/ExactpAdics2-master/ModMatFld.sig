174,0
T,ModMatFld_FldPadExact,ModMatFldElt_FldPadExact,1,StrPadExact
A,ModMatFld_FldPadExact,10,base_field,nrows,ncols,generators,degree,dimension,row_space,transpose_space,identity,zero
T,ModMatFldElt_FldPadExact,-,1,PadExactElt
A,ModMatFldElt_FldPadExact,7,negation,component,row,rows,eltseq,determinant,inverse
A,FldPadExact,1,matrix_space
S,KMatrixSpace,The space of `m x n` matrices over K,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,FldPadExact,,ModMatFld_FldPadExact,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,ModMatFld_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,ModMatFldElt_FldPadExact,,-38,-38,-38,-38,-38,-38
S,InterpolateEpochs,Interpolates between the given epochs,0,4,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,148,,0,0,ModMatFldElt_FldPadExact,,168,-38,-38,-38,-38,-38
S,Ngens,The number of generators of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,NumberOfNames,The number of generators of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Name,The ith generator of M,0,2,0,0,0,0,0,0,0,148,,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,.,The ith generator of M,0,2,0,0,0,0,0,0,0,148,,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,AssignNames,Assigns names to the generators of V,1,1,1,82,0,298,2,0,0,1,0,0,0,0,82,,1,1,ModMatFld_FldPadExact,,-38,-38,-38,-38,-38,-38
S,Matrix,The `nrows x ncols` matrix over K defined by cs,0,4,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,148,,0,0,FldPadExactElt,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Matrix,The `nrows x ncols` matrix with coefficients cs,1,2,1,82,0,FldPadExactElt,3,0,0,0,0,0,0,0,82,,0,0,148,,0,0,148,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Matrix,The matrix whose rows are given by each entry in cs,1,0,1,82,1,82,0,FldPadExactElt,1,0,0,0,0,0,0,0,82,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Matrix,The matrix whose rows are given by each entry in cs,1,0,1,82,0,ModTupFldElt_FldPadExact,1,0,0,0,0,0,0,0,82,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Zero,The zero matrix,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ZeroMatrix,The zero matrix,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ZeroMatrix,The zero matrix,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Identity,The identity matrix,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,IdentityMatrix,The identity matrix,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,IdentityMatrix,The identity matrix,0,2,0,0,0,0,0,0,0,148,,0,0,FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ScalarMatrix,The scalar matrix with x on the diagonal and zero elsewhere,0,2,0,0,0,0,0,0,0,-1,,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ScalarMatrix,The scalar matrix with x on the diagonal and zero elsewhere,0,3,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,ScalarMatrix,The scalar matrix with x on the diagonal and zero elsewhere,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,148,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,DiagonalMatrix,The diagonal matrix with diagonal entries given by cs,0,2,0,0,0,0,0,0,0,82,,0,0,ModMatFld_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,DiagonalMatrix,The diagonal matrix with diagonal entries given by cs,0,2,0,0,0,0,0,0,0,82,,0,0,FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,DiagonalMatrix,The diagonal matrix with diagonal entries given by cs,1,0,1,82,0,FldPadExactElt,1,0,0,0,0,0,0,0,82,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of M. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,-1,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of M. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of M. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,82,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of M. If so, also returns the coerced element",1,1,1,82,0,82,2,0,0,0,0,0,0,0,82,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of M. If so, also returns the coerced element",1,1,1,82,0,ModTupFldElt_FldPadExact,2,0,0,0,0,0,0,0,82,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,BaseField,The base field of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,BaseField,The base field of m,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,FldPadExact,-38,-38,-38,-38,-38
S,Nrows,The number of rows of elements of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Ncols,The number of columns of elements of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Degree,"If V is a subspace of `K^(m x n)`, returns `mn`. That is, the number of components in elements of M",0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Dimension,The dimension of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,148,-38,-38,-38,-38,-38
S,Generators,The generators of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,82,-38,-38,-38,-38,-38
S,ExistsCoveringStructure,True if there is a space containing both M1 and M2,0,2,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,0,0,ModMatFld_FldPadExact,,36,-1,-38,-38,-38,-38
S,CanChangeRing,True if the base ring of m can be changed to F,0,2,0,0,0,0,0,0,0,FldPadExact,,0,0,ModMatFldElt_FldPadExact,,36,ModMatFldElt_FldPadExact,-38,-38,-38,-38
S,ChangeRing,A copy of m with its base field changed to F,0,2,0,0,0,0,0,0,0,FldPadExact,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,RowSpace,The vector space of rows of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModTupFld_FldPadExact,-38,-38,-38,-38,-38
S,TransposeSpace,The space of transposes of elements of M,0,1,0,0,0,0,0,0,0,ModMatFld_FldPadExact,,ModMatFld_FldPadExact,-38,-38,-38,-38,-38
S,Nrows,Number of rows,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,148,-38,-38,-38,-38,-38
S,Ncols,Number of columns,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,148,-38,-38,-38,-38,-38
S,Eltseq,The components of m,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,Component,The jth component of the ith row of m,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,ModMatFldElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,@,The jth component of the ith row of m,0,3,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,148,,0,0,148,,FldPadExactElt,-38,-38,-38,-38,-38
S,Rows,The rows of m,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,82,-38,-38,-38,-38,-38
S,Row,The ith row of m,0,2,0,0,0,0,0,0,0,148,,0,0,ModMatFldElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,@,The ith row of m,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,148,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,-,Negation,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,-1,,0,0,ModMatFldElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,+,Addition,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,-1,,0,0,ModMatFldElt_FldPadExact,,-1,-38,-38,-38,-38,-38
S,-,Subtraction,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,&+,Sum,1,0,1,82,0,ModMatFldElt_FldPadExact,1,0,0,0,0,0,0,0,82,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Scalar multiplication,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Scalar multiplication,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,FldPadExactElt,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,/,Scalar division,0,2,0,0,0,0,0,0,0,FldPadExactElt,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Vector-matrix multiplication,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,ModTupFldElt_FldPadExact,,ModTupFldElt_FldPadExact,-38,-38,-38,-38,-38
S,*,Matrix multiplication,0,2,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,^,Matrix power. Negative powers are allowed for invertible matrices,0,2,0,0,0,0,0,0,0,148,,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Transpose,Transpose,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,ModMatFldElt_FldPadExact,-38,-38,-38,-38,-38
S,Determinant,Determinant,0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,FldPadExactElt,-38,-38,-38,-38,-38
S,IsDefinitelyInvertible,"True if m is definitely invertible. If so, also returns the inverse",0,1,0,0,0,0,0,0,0,ModMatFldElt_FldPadExact,,36,ModMatFldElt_FldPadExact,-38,-38,-38,-38
