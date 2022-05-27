174,0
T,PGGAlg_ResPriority,-,1,PGGAlg
T,PGGAlg_ResPriority_Null,-,1,PGGAlg_ResPriority
T,PGGAlg_ResPriority_Random,-,1,PGGAlg_ResPriority
T,PGGAlg_ResPriority_Reverse,-,1,PGGAlg_ResPriority
A,PGGAlg_ResPriority_Reverse,1,priority
T,PGGAlg_ResPriority_Key,-,1,PGGAlg_ResPriority
A,PGGAlg_ResPriority_Key,1,key
T,PGGAlg_ResPriority_Truncate,-,1,PGGAlg_ResPriority
A,PGGAlg_ResPriority_Truncate,2,priority,limit
S,PGGAlg_ResPriority_Null_Make,"The ""Null"" subgroups priority",0,0,0,0,0,0,0,PGGAlg_ResPriority_Null,-38,-38,-38,-38,-38
S,PGGAlg_ResPriority_Random_Make,"The ""Random"" subgroups priority",0,0,0,0,0,0,0,PGGAlg_ResPriority_Random,-38,-38,-38,-38,-38
S,PGGAlg_ResPriority_Reverse_Make,"The ""Reverse"" subgroups priority",0,0,0,0,0,0,0,PGGAlg_ResPriority_Reverse,-38,-38,-38,-38,-38
S,PGGAlg_ResPriority_Truncate_Make,"The ""truncate"" subgroups priority",0,1,0,0,0,0,0,0,0,148,,PGGAlg_ResPriority_Truncate,-38,-38,-38,-38,-38
S,PGGAlg_ResPriority_Key_Make,"The ""Key"" subgroups priority",0,1,0,0,0,0,0,0,0,PGGExpr,,PGGAlg_ResPriority_Key,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResPriority_Null,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResPriority_Random,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResPriority_Reverse,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResPriority_Truncate,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResPriority_Key,,-38,-38,-38,-38,-38,-38
S,Prioritize,Prioritizes X,0,3,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResPriority_Null,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,Prioritize,Prioritizes X,0,3,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResPriority_Reverse,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,Prioritize,Prioritizes X,0,3,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResPriority_Random,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,Prioritize,Prioritizes X,0,3,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResPriority_Key,,0,0,-1,,-1,-38,-38,-38,-38,-38
S,Prioritize,Prioritizes X,0,3,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResPriority_Truncate,,0,0,-1,,-1,-38,-38,-38,-38,-38
