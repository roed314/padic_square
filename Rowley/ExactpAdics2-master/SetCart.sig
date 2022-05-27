174,0
T,SetCart_PadExact,Tup_PadExact,1,StrPadExact
A,SetCart_PadExact,1,components
T,Tup_PadExact,-,1,PadExactElt
A,Tup_PadExact,2,to_tuple,to_sequence
S,ExactpAdics_CartesianProduct,The cartesian product of the components,0,1,0,0,0,0,0,0,0,303,,SetCart_PadExact,-38,-38,-38,-38,-38
S,CartesianPower,The nth cartesian power of S,0,2,0,0,0,0,0,0,0,148,,0,0,StrPadExact,,SetCart_PadExact,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,SetCart_PadExact,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,2,0,0,1,0,0,0,0,298,,0,0,Tup_PadExact,,-38,-38,-38,-38,-38,-38
S,InterpolateEpochs,Interpolates between the given epochs,0,4,0,0,0,0,0,0,0,303,,0,0,148,,0,0,148,,0,0,Tup_PadExact,,168,-38,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of C. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,-1,,0,0,SetCart_PadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of C. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,Tup_PadExact,,0,0,SetCart_PadExact,,36,-1,-38,-38,-38,-38
S,IsCoercible,"True if X is coercible to an element of C. If so, also returns the coerced element",0,2,0,0,0,0,0,0,0,303,,0,0,SetCart_PadExact,,36,-1,-38,-38,-38,-38
S,NumberOfComponents,The number of components of C,0,1,0,0,0,0,0,0,0,SetCart_PadExact,,148,-38,-38,-38,-38,-38
S,@,The ith component of C,0,2,0,0,0,0,0,0,0,SetCart_PadExact,,0,0,148,,StrPadExact,-38,-38,-38,-38,-38
S,Components,The components of C,0,1,0,0,0,0,0,0,0,SetCart_PadExact,,303,-38,-38,-38,-38,-38
S,#,The number of elements of T,0,1,0,0,0,0,0,0,0,Tup_PadExact,,148,-38,-38,-38,-38,-38
S,@,The ith element of T,0,2,0,0,0,0,0,0,0,Tup_PadExact,,0,0,148,,PadExactElt,-38,-38,-38,-38,-38
S,ToTuple,Converts T to a regular tuple,0,1,0,0,0,0,0,0,0,Tup_PadExact,,303,-38,-38,-38,-38,-38
S,ToSequence,Converts T to a sequence,0,1,0,0,0,0,0,0,0,Tup_PadExact,,82,-38,-38,-38,-38,-38
