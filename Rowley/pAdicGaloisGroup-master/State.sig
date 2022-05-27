174,1
T,PGGState,-,0
A,PGGState,3,timer,factorization_func,roots_func
V,PGG_GlobalTimer,1
S,PGG_Factorization,"Factorization of f as a sequence of factors (assumes f is squarefree). Also returns a sequence of certificates. If Extensions is true, then the certificates have the `Extension` field assigned. If Lift is false, the factors are not Hensel lifted, and so may have low precision",1,0,1,312,0,400,1,0,0,0,0,0,0,0,312,,82,82,-38,-38,-38,-38
S,PGG_Roots,"Roots of f as a sequence (assumes f is squarefree). If Lift is false, the roots are not Hensel lifted, and so may have low precision",1,0,1,312,0,400,1,0,0,0,0,0,0,0,312,,82,-38,-38,-38,-38,-38
S,PGG_UseBuiltinRoots,Uses the builtin roots algorithm,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_UseExactpAdicsRoots,Uses the new roots algorithm. Requires the ExactpAdics package to be attached,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_UseBuiltinFactorization,Uses the builtin factorization algorithm,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_UseExactpAdicsFactorization,Uses the new factorization algorithm. Requires the ExactpAdics package to be attached,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_HasGlobalTimer,True if there is a global timer set,0,0,0,0,0,0,0,36,PGGTimer,-38,-38,-38,-38
S,PGG_GlobalTimer,The global timer,0,0,0,0,0,0,0,PGGTimer,-38,-38,-38,-38,-38
S,PGG_StopGlobalTimer,Stops the global timer,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_StartGlobalTimer,Starts the global timer,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Push,Push,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_GetLabel,Gets the current label,0,0,0,0,0,0,0,82,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_SetLabel,Sets the current label,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Pop,Pop,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Pop,Pop,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Pop,Pop,0,2,0,0,1,0,0,0,0,-1,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Swap,Swap,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_Log,Log,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,PGG_GlobalTimer_PrintTree,PrintTree,0,0,0,0,1,0,0,-38,-38,-38,-38,-38,-38
