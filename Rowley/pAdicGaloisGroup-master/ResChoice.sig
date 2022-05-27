174,0
T,PGGAlg_ResChoice,-,1,PGGAlg
T,PGGAlg_ResChoice_Tranche,-,1,PGGAlg_ResChoice
A,PGGAlg_ResChoice_Tranche,2,tranche_alg,priority
T,PGGAlg_ResChoice_Stream,-,1,PGGAlg_ResChoice
A,PGGAlg_ResChoice_Stream,2,stream_alg,limit
T,PGGAlgState_ResChoice,-,1,PGGAlgState
A,PGGAlgState_ResChoice,1,parent
T,PGGAlgState_ResChoice_Tranche,-,1,PGGAlgState_ResChoice
A,PGGAlgState_ResChoice_Tranche,2,tranche_state,overgroup
T,PGGAlgState_ResChoice_Stream,-,1,PGGAlgState_ResChoice
A,PGGAlgState_ResChoice_Stream,2,stream_state,overgroup
S,PGGAlg_ResChoice_Tranche_Make,"The ""Tranche"" subgroups choice",0,0,0,0,0,0,0,PGGAlg_ResChoice_Tranche,-38,-38,-38,-38,-38
S,PGGAlg_ResChoice_Stream_Make,"The ""Stream"" subgroups choice",0,1,0,0,0,0,0,0,0,PGGAlg_Stream,,PGGAlg_ResChoice_Stream,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResChoice_Tranche,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResChoice_Stream,,-38,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResChoice_Tranche,,PGGAlgState_ResChoice_Tranche,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResChoice_Stream,,PGGAlg_ResChoice_Stream,-38,-38,-38,-38,-38
S,Reset,Starts the choice from the beginning again,0,1,0,0,1,0,0,0,0,PGGAlgState_ResChoice_Tranche,,-38,-38,-38,-38,-38,-38
S,Reset,Starts the choice from the beginning again,0,1,0,0,1,0,0,0,0,PGGAlgState_ResChoice_Stream,,-38,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`parent)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResChoice,,0,0,PGGAlg_ResChoice,,PGGAlgState_ResChoice,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`parent)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResChoice_Tranche,,0,0,PGGAlg_ResChoice_Tranche,,PGGAlgState_ResChoice_Tranche,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`parent)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResChoice_Stream,,0,0,PGGAlg_ResChoice_Stream,,PGGAlgState_ResChoice_Stream,-38,-38,-38,-38,-38
S,HasSubgroup,Selects a subgroup,0,1,0,0,0,0,0,0,0,PGGAlgState_ResChoice_Tranche,,36,224,-38,-38,-38,-38
S,HasSubgroup,Selects a subgroup,0,1,0,0,0,0,0,0,0,PGGAlgState_ResChoice_Stream,,36,224,-38,-38,-38,-38
S,Subgroup,The next subgroup,0,1,0,0,0,0,0,0,0,PGGAlgState_ResChoice,,224,-38,-38,-38,-38,-38
S,MarkSpecial,"Marks the current state of s as ""special""",0,1,0,0,1,0,0,0,0,PGGAlgState_ResChoice_Tranche,,-38,-38,-38,-38,-38,-38
