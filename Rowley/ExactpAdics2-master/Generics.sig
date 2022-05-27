174,0
T,AnyPadExact,-,0
A,AnyPadExact,7,id,approximations,dependencies,get_approximation,min_epoch,max_epoch,internal_data
T,StrPadExact,PadExactElt,1,AnyPadExact
T,PadExactElt,-,1,AnyPadExact
A,PadExactElt,1,parent
T,ValPadExact,-,0
S,BringToEpoch,Brings x to epoch n,0,2,0,0,1,0,0,0,0,148,,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,BringToEpoch,Brings each `x` in `xs` to epoch `n`,1,0,1,82,0,AnyPadExact,2,0,0,1,0,0,0,0,148,,0,0,82,,-38,-38,-38,-38,-38,-38
S,CanBringToEpoch,"Tries to bring x to epoch n. Returns `true` if successful, or `false` if the max_epoch of something is exceeded",0,2,0,0,0,0,0,0,0,148,,0,0,AnyPadExact,,36,-38,-38,-38,-38,-38
S,CanBringToEpoch,"Tries to bring each `x` in `xs` to epoch `n`. Returns `true` if successful, or `false` if the max_epoch of something is exceeded",1,0,1,82,0,AnyPadExact,2,0,0,0,0,0,0,0,148,,0,0,82,,36,-38,-38,-38,-38,-38
S,EpochApproximation,Returns sequence of approximations of xs at epoch n,1,0,1,82,0,AnyPadExact,2,0,0,0,0,0,0,0,148,,0,0,82,,82,-38,-38,-38,-38,-38
S,EpochApproximation,The approximation of x at epoch n,0,2,0,0,0,0,0,0,0,148,,0,0,AnyPadExact,,-1,-38,-38,-38,-38,-38
S,Init_AnyPadExact,Initializes x,0,1,0,0,1,0,0,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,Init_PadExactElt,Initializes x,0,1,0,0,1,0,0,0,0,PadExactElt,,-38,-38,-38,-38,-38,-38
S,Init,Initializes x,0,1,0,0,1,0,0,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,Init,Initializes x,0,1,0,0,1,0,0,0,0,PadExactElt,,-38,-38,-38,-38,-38,-38
S,Parent,The parent of x,0,1,0,0,0,0,0,0,0,PadExactElt,,StrPadExact,-38,-38,-38,-38,-38
S,BestApproximation,The best approximation to x,0,1,0,0,0,0,0,0,0,AnyPadExact,,-1,-38,-38,-38,-38,-38
S,eq,Equality,0,2,0,0,0,0,0,0,0,StrPadExact,,0,0,StrPadExact,,148,-38,-38,-38,-38,-38
S,ExistsEpochWith,"True if there is an epoch of x such that `pred(x)` is true. If so, returns it",0,2,0,0,0,0,0,0,0,-1,,0,0,AnyPadExact,,36,148,-38,-38,-38,-38
S,ExistsEpochWithApproximation,"True if there is an epoch of x such that `pred(xx)` is true, where `xx` is the corresponding approximation. If so, returns it. False if `max_epoch` is ever exceeded",0,2,0,0,0,0,0,0,0,-1,,0,0,AnyPadExact,,36,148,-38,-38,-38,-38
S,FirstEpochWithApproximation,The first epoch whose approximation satisfies the predicate,0,2,0,0,0,0,0,0,0,-1,,0,0,AnyPadExact,,148,-38,-38,-38,-38,-38
S,IncreaseEpochUntil,Increases the epoch of x until `pred(x)` is true,0,2,0,0,1,0,0,0,0,-1,,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,IncreaseAbsolutePrecision,Increases the absolute precision of x to at least apr,0,2,0,0,1,0,0,0,0,-1,,0,0,PadExactElt,,-38,-38,-38,-38,-38,-38
S,Approximation,An approximation of x with absolute precision at least apr,0,2,0,0,0,0,0,0,0,-1,,0,0,PadExactElt,,-1,-38,-38,-38,-38,-38
S,InterpolateUpTo,Replaces the approximations of x below n with their interpolations,0,2,0,0,1,0,0,0,0,148,,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,InterpolateUpToFirstEpochWithApproximation,"Replaces the approximations of x below n with their interpolations, where n is the first epoch whose approximation satisfies the predicate. Returns the new and old values of n",0,2,0,0,0,0,0,0,0,-1,,0,0,AnyPadExact,,148,148,-38,-38,-38,-38
S,SetApproximation,Sets the approximation of x at epoch n to xx,0,3,0,0,1,0,0,0,0,-1,,0,0,148,,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
V,ExactpAdics_WithDependencies,1
S,WithDependencies,A copy of x with the direct dependencies ds,0,2,0,0,0,0,0,0,0,168,,0,0,AnyPadExact,,AnyPadExact,-38,-38,-38,-38,-38
S,InterpolateEpochs,The approximations for epochs [n1+1..n2-1],0,4,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,148,,0,0,AnyPadExact,,168,-38,-38,-38,-38,-38
S,GetApproximation,Gets an approximation using get,0,3,0,0,0,0,0,0,0,168,,0,0,148,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,GetApproximation,Gets an approximation using get,0,3,0,0,0,0,0,0,0,168,,0,0,148,,0,0,41,,-1,-38,-38,-38,-38,-38
S,SetApproximationHook,Called by `SetApproximation`,0,3,0,0,1,0,0,1,1,-38,,0,0,148,,0,0,AnyPadExact,,-38,-38,-38,-38,-38,-38
S,IsValidApproximation,"True if xx is a valid approximation of x at epoch n. If not, also returns an error message",0,3,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,AnyPadExact,,36,298,-38,-38,-38,-38
S,IsValidApproximation,"True if xx is a valid approximation of x at epoch n. If not, also returns an error message",0,3,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,PadExactElt,,36,298,-38,-38,-38,-38
S,UninitializedCopy,An uninitialized copy of x,0,1,0,0,0,0,0,0,0,AnyPadExact,,AnyPadExact,-38,-38,-38,-38,-38
S,UninitializedCopy,An uninitialized copy of x,0,1,0,0,0,0,0,0,0,StrPadExact,,StrPadExact,-38,-38,-38,-38,-38
S,UninitializedCopy,An uninitialized copy of x,0,1,0,0,0,0,0,0,0,PadExactElt,,PadExactElt,-38,-38,-38,-38,-38
