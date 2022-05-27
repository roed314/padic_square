174,0
T,PGGAlg_ResEval,-,1,PGGAlg
T,PGGAlg_ResEval_Global,-,1,PGGAlg_ResEval
A,PGGAlg_ResEval_Global,1,model
T,PGGAlgState_ResEval_Global,-,1,PGGAlgState
A,PGGAlgState_ResEval_Global,6,base_field_model,pol_model,precision,base_complex_embeddings,complex_roots,overgroup_embedding
T,PGGAlg_ResEval_Global_Model,-,1,PGGAlg
T,PGGAlg_ResEval_Global_Model_Symmetric,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_Symmetric,1,galois_group_alg
T,PGGAlg_ResEval_Global_Model_Factors,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_Factors,1,next
T,PGGAlg_ResEval_Global_Model_Tower,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_Tower,1,next
T,PGGAlg_ResEval_Global_Model_RamTower,-,1,PGGAlg_ResEval_Global_Model_Tower
T,PGGAlg_ResEval_Global_Model_D4Tower,-,1,PGGAlg_ResEval_Global_Model_Tower
T,PGGAlg_ResEval_Global_Model_Select,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_Select,2,predicates,models
T,PGGAlg_ResEval_Global_Model_RootOfUnity,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_RootOfUnity,2,minimize,complement
T,PGGAlg_ResEval_Global_Model_RootOfUniformizer,-,1,PGGAlg_ResEval_Global_Model
T,PGGAlg_ResEval_Global_Model_PthRoots,-,1,PGGAlg_ResEval_Global_Model
T,PGGAlg_ResEval_Global_Model_Cheat,-,1,PGGAlg_ResEval_Global_Model
A,PGGAlg_ResEval_Global_Model_Cheat,1,next
S,Start,Starts the algorithm and returns its state,0,2,0,0,0,0,0,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global,,PGGAlgState_ResEval_Global,-38,-38,-38,-38,-38
S,OvergroupEmbedding,The overgroup we can evaluate resolvents in,0,1,0,0,0,0,0,0,0,PGGAlgState_ResEval_Global,,PGGHomGrpPerm,-38,-38,-38,-38,-38
S,ComplexRoots,"A sequence of complex roots, one for each embedding of the base field. Also returns the base field embeddings",0,2,0,0,0,0,0,0,0,148,,0,0,PGGAlgState_ResEval_Global,,82,82,-38,-38,-38,-38
S,Resolvent,The resolvent associated to the invariant I,0,3,0,0,0,0,0,0,0,224,,0,0,-1,,0,0,PGGAlgState_ResEval_Global,,PGGPol,-38,-38,-38,-38,-38
S,Resolvent,The resolvent associated to the invariant I,0,4,0,0,0,0,0,0,0,224,,0,0,-1,,0,0,PGGGloMod,,0,0,PGGAlgState_ResEval_Global,,PGGPol,-38,-38,-38,-38,-38
S,Resolvent,The resolvent associated to the invariant I,0,4,0,0,0,0,0,0,0,224,,0,0,-1,,0,0,PGGGloMod_UPol_Cheat,,0,0,PGGAlgState_ResEval_Global,,PGGPol,-38,-38,-38,-38,-38
S,Resolvent,The resolvent associated to the invariant I,0,5,0,0,0,0,0,0,0,224,,0,0,-1,,0,0,PGGPolGrp,,0,0,PGGGloMod_UPol_Cheat,,0,0,PGGAlgState_ResEval_Global,,PGGPol,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_Symmetric_Make,"The ""Symmetric"" conjugacy algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_Symmetric,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_Factors_Make,"The ""Factors"" conjugacy algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_Factors,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_RamTower_Make,"The ""RamTower"" conjugacy algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_RamTower,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_D4Tower_Make,"The ""D4Tower"" conjugacy algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_D4Tower,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_Select_Make,"The ""Select"" global model algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_Select,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_RootOfUnity_Make,"The ""RootOfUnity"" global model algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_RootOfUnity,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_RootOfUniformizer_Make,"The ""RootOfUniformizer"" global model algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_RootOfUniformizer,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_PthRoots_Make,"The ""PthRoots"" global model algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_PthRoots,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Model_Cheat_Make,"The ""Cheat"" global model algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global_Model_Cheat,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_Symmetric,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_Factors,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_RamTower,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_D4Tower,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_Select,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_RootOfUnity,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_RootOfUniformizer,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_PthRoots,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global_Model_Cheat,,-38,-38,-38,-38,-38,-38
S,PGGAlg_ResEval_Global_Make,"The ""Absolute"" resolvent evaluation algorithm",0,0,0,0,0,0,0,PGGAlg_ResEval_Global,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResEval_Global,,-38,-38,-38,-38,-38,-38
