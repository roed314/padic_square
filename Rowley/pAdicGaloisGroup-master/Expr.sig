174,0
T,PGGExpr,-,0
A,PGGExpr,1,free_vars
T,PGGExpr_FreeVar,-,1,PGGExpr
A,PGGExpr_FreeVar,1,name
T,PGGExpr_Unop,-,1,PGGExpr
A,PGGExpr_Unop,5,op,arg1,name,prefix,format
T,PGGExpr_Binop,-,1,PGGExpr
A,PGGExpr_Binop,6,op,arg1,arg2,name,infix,format
T,PGGExpr_Seqop,-,1,PGGExpr
A,PGGExpr_Seqop,4,op,args,name,infix
T,PGGExpr_Const,-,1,PGGExpr
A,PGGExpr_Const,1,val
T,PGGExpr_All,-,1,PGGExpr
A,PGGExpr_All,1,args
T,PGGExpr_Any,-,1,PGGExpr
A,PGGExpr_Any,1,args
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_FreeVar,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_Unop,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_Binop,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_Seqop,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_Const,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_All,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGExpr_Any,,-38,-38,-38,-38,-38,-38
S,FreeVariables,The free variables in the expression,0,1,0,0,0,0,0,0,0,PGGExpr,,83,-38,-38,-38,-38,-38
S,PGG_Expression_FreeVariable,The free variable with the given name,0,1,0,0,0,0,0,0,0,298,,PGGExpr_FreeVar,-38,-38,-38,-38,-38
S,PGG_Expression_UnOp,The unary operation op (a function or intrinsic) with argument arg1,0,2,0,0,0,0,0,0,0,PGGExpr,,0,0,-1,,PGGExpr_Unop,-38,-38,-38,-38,-38
S,PGG_Expression_BinOp,The binary operation op (a function or intrinsic) with arguments arg1 and arg2,0,3,0,0,0,0,0,0,0,PGGExpr,,0,0,PGGExpr,,0,0,-1,,PGGExpr_Binop,-38,-38,-38,-38,-38
S,PGG_Expression_SeqOp,"The operation op with one argument, the sequence args",0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,PGGExpr_Seqop,-38,-38,-38,-38,-38
S,PGG_Expression_Const,The constant expression val,0,1,0,0,0,0,0,0,0,-1,,PGGExpr_Const,-38,-38,-38,-38,-38
S,PGG_Expression_All,The expression evaluating to true if all its arguments are true,0,1,0,0,0,0,0,0,0,-1,,PGGExpr_All,-38,-38,-38,-38,-38
S,PGG_Expression_Any,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,-1,,PGGExpr_Any,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_FreeVar,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_Unop,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_Binop,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_Const,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_Seqop,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_All,,83,-38,-38,-38,-38,-38
S,_FreeVariables,The expression evaluating to true if any of its arguments are true,0,1,0,0,0,0,0,0,0,PGGExpr_Any,,83,-38,-38,-38,-38,-38
S,Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,41,,0,0,PGGExpr,,-1,-38,-38,-38,-38,-38
S,Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,457,,0,0,PGGExpr,,-1,-38,-38,-38,-38,-38
S,EvaluateLazy,Evaluates x,0,2,0,0,0,0,0,0,0,457,,0,0,PGGExpr,,-1,-1,-38,-38,-38,-38
S,EvaluateLazy,Evaluates x,0,2,0,0,0,0,0,0,0,41,,0,0,PGGExpr,,-1,-1,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_FreeVar,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_Unop,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_Binop,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_Const,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_Seqop,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_All,,-1,-38,-38,-38,-38,-38
S,_Evaluate,Evaluates x,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGExpr_Any,,-1,-38,-38,-38,-38,-38
