174,0
T,PGGAlg_Tranche,-,1,PGGAlg
A,PGGAlg_Tranche,1,verbosity
T,PGGAlg_Tranche_All,-,1,PGGAlg_Tranche
T,PGGAlg_Tranche_Index,-,1,PGGAlg_Tranche
A,PGGAlg_Tranche_Index,6,filter,sort_key,static_filter,dynamic_filter,take,dedupe
T,PGGAlg_Tranche_OrbitIndex,-,1,PGGAlg_Tranche
A,PGGAlg_Tranche_OrbitIndex,6,filter,sort_key,static_filter,dynamic_filter,take,dedupe
T,PGGAlg_Tranche_Modify,-,1,PGGAlg_Tranche
A,PGGAlg_Tranche_Modify,1,inner
T,PGGAlg_Tranche_Tuples,-,1,PGGAlg_Tranche_Modify
A,PGGAlg_Tranche_Tuples,2,length,random_limit
T,PGGAlg_Tranche_Shuffle,-,1,PGGAlg_Tranche_Modify
T,PGGAlg_Tranche_Truncate,-,1,PGGAlg_Tranche_Modify
A,PGGAlg_Tranche_Truncate,1,length
T,PGGAlgState_Tranche,-,1,PGGAlgState
A,PGGAlgState_Tranche,2,overgroup,special_tranches
T,PGGAlgState_Tranche_All,-,1,PGGAlgState_Tranche
A,PGGAlgState_Tranche_All,2,all_groups,tranche
T,PGGAlgState_Tranche_Index,-,1,PGGAlgState_Tranche
A,PGGAlgState_Tranche_Index,6,all_indices,indices,ii,all_groups,tranche,subgroup_classes
T,PGGAlgState_Tranche_OrbitIndex,-,1,PGGAlgState_Tranche
A,PGGAlgState_Tranche_OrbitIndex,9,all_indices,indices,ii,si,all_groups,all_stabilizers,cur_stabilizers,tranche,subgroup_classes
T,PGGAlgState_Tranche_Modify,-,1,PGGAlgState_Tranche
A,PGGAlgState_Tranche_Modify,2,inner_state,tranche
T,PGGAlg_TrancheTake,-,1,PGGAlg
T,PGGAlg_TrancheTake_All,-,1,PGGAlg_TrancheTake
T,PGGAlg_TrancheTake_FrattiniBasis,-,1,PGGAlg_TrancheTake
T,PGGAlg_TrancheTake_Random,-,1,PGGAlg_TrancheTake
A,PGGAlg_TrancheTake_Random,3,limit,new_tries,random_tries
T,PGGAlgState_TrancheItem,-,0
A,PGGAlgState_TrancheItem,4,state,idx,tidx,subgroup
S,PGGAlg_Tranche_All_Make,"The ""All"" subgroups choice",0,0,0,0,0,0,0,PGGAlg_Tranche_All,-38,-38,-38,-38,-38
S,PGGAlg_Tranche_Index_Make,"The ""Index"" subgroups choice",0,0,0,0,0,0,0,PGGAlg_Tranche_Index,-38,-38,-38,-38,-38
S,PGGAlg_Tranche_OrbitIndex_Make,"The ""Index"" subgroups choice",0,0,0,0,0,0,0,PGGAlg_Tranche_OrbitIndex,-38,-38,-38,-38,-38
S,PGGAlg_Tranche_Tuples_Make,"The ""Tuples"" subgroups choice",0,2,0,0,0,0,0,0,0,PGGAlg_Tranche,,0,0,148,,PGGAlg_Tranche_Tuples,-38,-38,-38,-38,-38
S,PGGAlg_Tranche_Shuffle_Make,"The ""Shuffle"" subgroups choice",0,1,0,0,0,0,0,0,0,PGGAlg_Tranche,,PGGAlg_Tranche_Shuffle,-38,-38,-38,-38,-38
S,PGGAlg_Tranche_Truncate_Make,"The ""Truncate"" subgroups choice",0,2,0,0,0,0,0,0,0,PGGAlg_Tranche,,0,0,148,,PGGAlg_Tranche_Truncate,-38,-38,-38,-38,-38
S,PGGAlg_TrancheTake_All_Make,Takes all subgroups,0,0,0,0,0,0,0,PGGAlg_TrancheTake_All,-38,-38,-38,-38,-38
S,PGGAlg_TrancheTake_FrattiniBasis_Make,Takes just the subgroups generated from a basis of the Frattini quotient,0,0,0,0,0,0,0,PGGAlg_TrancheTake_FrattiniBasis,-38,-38,-38,-38,-38
S,PGGAlg_TrancheTake_Random_Make,Takes some randomly generated groups,0,1,0,0,0,0,0,0,0,148,,PGGAlg_TrancheTake_Random,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_All,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_Index,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_OrbitIndex,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_Tuples,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_Shuffle,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Tranche_Truncate,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_TrancheTake_All,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_TrancheTake_FrattiniBasis,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_TrancheTake_Random,,-38,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGGrpPerm,,0,0,PGGAlg_Tranche_All,,PGGAlgState_Tranche_All,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGGrpPerm,,0,0,PGGAlg_Tranche_Index,,PGGAlgState_Tranche_Index,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGGrpPerm,,0,0,PGGAlg_Tranche_OrbitIndex,,PGGAlgState_Tranche_OrbitIndex,-38,-38,-38,-38,-38
S,Start,Starts alg,0,2,0,0,0,0,0,0,0,PGGGrpPerm,,0,0,PGGAlg_Tranche_Modify,,PGGAlgState_Tranche_Modify,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche,,0,0,PGGAlg_Tranche,,PGGAlgState_Tranche,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche_All,,0,0,PGGAlg_Tranche_All,,PGGAlgState_Tranche_All,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche_Index,,0,0,PGGAlg_Tranche_Index,,PGGAlgState_Tranche_Index,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche_OrbitIndex,,0,0,PGGAlg_Tranche_OrbitIndex,,PGGAlgState_Tranche_OrbitIndex,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,0,0,PGGAlg_Tranche_Modify,,PGGAlgState_Tranche_Modify,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche,,0,0,PGGAlg_Tranche_Modify,,PGGAlgState_Tranche_Modify,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,0,0,PGGAlg_Tranche,,PGGAlgState_Tranche,-38,-38,-38,-38,-38
S,Reset,Starts choosing from the beginning again,0,1,0,0,1,0,0,0,0,PGGAlgState_Tranche,,-38,-38,-38,-38,-38,-38
S,TrancheIndex,The index of the current tranche,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_All,,-1,-38,-38,-38,-38,-38
S,TrancheIndex,The index of the current tranche,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Index,,-1,-38,-38,-38,-38,-38
S,TrancheIndex,The index of the current tranche,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_OrbitIndex,,-1,-38,-38,-38,-38,-38
S,TrancheIndex,The index of the current tranche,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,-1,-38,-38,-38,-38,-38
S,InitialState,Gets the initial state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_All,,-1,-38,-38,-38,-38,-38
S,InitialState,Gets the initial state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Index,,-1,-38,-38,-38,-38,-38
S,InitialState,Gets the initial state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_OrbitIndex,,-1,-38,-38,-38,-38,-38
S,InitialState,Gets the initial state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,-1,-38,-38,-38,-38,-38
S,CurrentState,Gets the current state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_All,,-1,-38,-38,-38,-38,-38
S,CurrentState,Gets the current state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Index,,-1,-38,-38,-38,-38,-38
S,CurrentState,Gets the current state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_OrbitIndex,,-1,-38,-38,-38,-38,-38
S,CurrentState,Gets the current state of the choice,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,-1,-38,-38,-38,-38,-38
S,SetState,Sets the state of the choice,0,2,0,0,1,0,0,0,0,-1,,0,0,PGGAlgState_Tranche_All,,-38,-38,-38,-38,-38,-38
S,SetState,Sets the state of the choice,0,2,0,0,1,0,0,0,0,-1,,0,0,PGGAlgState_Tranche_Index,,-38,-38,-38,-38,-38,-38
S,SetState,Sets the state of the choice,0,2,0,0,1,0,0,0,0,-1,,0,0,PGGAlgState_Tranche_OrbitIndex,,-38,-38,-38,-38,-38,-38
S,SetState,Sets the state of the choice,0,2,0,0,1,0,0,0,0,-1,,0,0,PGGAlgState_Tranche_Modify,,-38,-38,-38,-38,-38,-38
S,Tranche,Gets a tranche of potential subgroups,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche,,82,-38,-38,-38,-38,-38
S,HasTranche,Gets a tranche of potential subgroups,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_All,,36,82,-38,-38,-38,-38
S,HasTranche,Gets a tranche of potential subgroups,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Index,,36,82,-38,-38,-38,-38
S,HasTranche,Gets a tranche of potential subgroups,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_OrbitIndex,,36,82,-38,-38,-38,-38
S,HasTranche,Gets a tranche of potential subgroups,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche_Modify,,36,82,-38,-38,-38,-38
S,Modify,"Givens a sequence of subgroups, returns the modified sequence of subgroups",0,2,0,0,0,0,0,0,0,82,,0,0,PGGAlg_Tranche_Modify,,82,-38,-38,-38,-38,-38
S,Modify,"Givens a sequence of subgroups, returns the modified sequence of subgroups",0,2,0,0,0,0,0,0,0,82,,0,0,PGGAlg_Tranche_Tuples,,82,-38,-38,-38,-38,-38
S,Modify,"Givens a sequence of subgroups, returns the modified sequence of subgroups",0,2,0,0,0,0,0,0,0,82,,0,0,PGGAlg_Tranche_Shuffle,,82,-38,-38,-38,-38,-38
S,Modify,"Givens a sequence of subgroups, returns the modified sequence of subgroups",0,2,0,0,0,0,0,0,0,82,,0,0,PGGAlg_Tranche_Truncate,,82,-38,-38,-38,-38,-38
S,SubgroupsOfIndex,Subgroups of G of index n,0,4,0,0,0,0,0,0,0,148,,0,0,224,,0,0,_PGGSetSubgrpcls,,0,0,PGGAlg_TrancheTake,,82,-38,-38,-38,-38,-38
S,SubgroupsOfIndex,Subgroups of G of index n,0,4,0,0,0,0,0,0,0,148,,0,0,224,,0,0,_PGGSetSubgrpcls,,0,0,PGGAlg_TrancheTake_All,,82,-38,-38,-38,-38,-38
S,SubgroupsOfIndex,Subgroups of G of index n,0,4,0,0,0,0,0,0,0,148,,0,0,224,,0,0,_PGGSetSubgrpcls,,0,0,PGGAlg_TrancheTake_FrattiniBasis,,82,-38,-38,-38,-38,-38
S,SubgroupsOfIndex,Subgroups of G of index n,0,4,0,0,0,0,0,0,0,148,,0,0,224,,0,0,_PGGSetSubgrpcls,,0,0,PGGAlg_TrancheTake_Random,,82,-38,-38,-38,-38,-38
S,MakeTranche,Makes a new tranche,0,2,0,0,0,0,0,0,0,82,,0,0,PGGAlgState_Tranche,,82,-38,-38,-38,-38,-38
S,TrancheItem,An item in a tranche,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,PGGAlgState_Tranche,,PGGAlgState_TrancheItem,-38,-38,-38,-38,-38
S,CosetAction,The coset action of i on G,0,2,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,0,0,224,,-1,-38,-38,-38,-38,-38
S,PGG_CosetAction,The coset action of i on G,0,2,0,0,0,0,0,0,0,-1,,0,0,224,,-1,-38,-38,-38,-38,-38
S,PGG_CosetAction,The coset action of i on G,0,2,0,0,0,0,0,0,0,224,,0,0,224,,-1,-38,-38,-38,-38,-38
S,PGG_CosetAction,The coset action of i on G,1,1,1,82,0,224,2,0,0,0,0,0,0,0,82,,0,0,224,,-1,-38,-38,-38,-38,-38
S,Id,Identifies this tranche item uniquely,0,1,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,-1,-38,-38,-38,-38,-38
S,Hash,Hash,0,1,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,-1,-38,-38,-38,-38,-38
S,eq,Equality,0,2,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,0,0,PGGAlgState_TrancheItem,,-1,-38,-38,-38,-38,-38
S,ForgetTrancheItem,Don't consider this item again,0,2,0,0,1,0,0,0,0,148,,0,0,PGGAlgState_Tranche,,-38,-38,-38,-38,-38,-38
S,ForgetTranche,Don't consider this item again,0,1,0,0,1,0,0,0,0,PGGAlgState_Tranche,,-38,-38,-38,-38,-38,-38
S,Group,The group,0,1,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,224,-38,-38,-38,-38,-38
S,Forget,Don't consider this item again,0,1,0,0,1,0,0,0,0,PGGAlgState_TrancheItem,,-38,-38,-38,-38,-38,-38
S,MarkSpecialTranche,"Mark this tranche as ""special""",0,1,0,0,1,0,0,0,0,PGGAlgState_TrancheItem,,-38,-38,-38,-38,-38,-38
S,MarkSpecialTranche,"Marks the tranche tidx of s as ""special""",0,2,0,0,1,0,0,0,0,-1,,0,0,PGGAlgState_Tranche,,-38,-38,-38,-38,-38,-38
S,MarkSpecialTranche,"Marks the current tranche of s as ""special""",0,1,0,0,1,0,0,0,0,PGGAlgState_Tranche,,-38,-38,-38,-38,-38,-38
S,SpecialTranches,The special tranches of s,0,1,0,0,0,0,0,0,0,PGGAlgState_Tranche,,151,-38,-38,-38,-38,-38
