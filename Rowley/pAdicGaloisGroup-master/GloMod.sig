174,0
T,PGGGloMod,-,0
T,PGGGloMod_Fld,-,1,PGGGloMod
A,PGGGloMod_Fld,5,local_field,global_field,embeddings,approximations,primes
T,PGGGloMod_UPol,-,1,PGGGloMod
A,PGGGloMod_UPol,2,local_pol,global_pol
T,PGGGloMod_Ext,-,1,PGGGloMod
A,PGGGloMod_Ext,1,base_model
T,PGGGloMod_FldExt,-,2,PGGGloMod_Fld,PGGGloMod_Ext
T,PGGGloMod_Rational,-,1,PGGGloMod_Fld
T,PGGGloMod_Symmetric,-,2,PGGGloMod_FldExt,PGGGloMod_UPol
T,PGGGloMod_Factors,-,1,PGGGloMod_Ext
A,PGGGloMod_Factors,1,factors
T,PGGGloMod_Tower,-,1,PGGGloMod_FldExt
A,PGGGloMod_Tower,3,tower,global_pol,global_pol_root
T,PGGGloMod_RootOfUnity,-,2,PGGGloMod_FldExt,PGGGloMod_UPol
A,PGGGloMod_RootOfUnity,4,primitive_root,aut_powers,sum_powers,zeta_pol
T,PGGGloMod_RootOfUniformizer,-,2,PGGGloMod_FldExt,PGGGloMod_UPol
A,PGGGloMod_RootOfUniformizer,3,degree,global_uniformizer,global_zeta_pol
T,PGGGloMod_PthRoots,-,2,PGGGloMod_FldExt,PGGGloMod_UPol
A,PGGGloMod_PthRoots,6,p,n,local_gens,global_gens,global_root_coeffs,global_zetap
T,PGGGloMod_UPol_Cheat,-,1,PGGGloMod_UPol
A,PGGGloMod_UPol_Cheat,2,model,overgroup_embedding
S,Approximation,Takes a global approximation to x across all embeddings,0,3,0,0,0,0,0,0,0,148,,0,0,PGGFldElt,,0,0,PGGGloMod_Fld,,-25,-38,-38,-38,-38,-38
S,Approximation,Takes a global approximation to x across all embeddings,0,3,0,0,0,0,0,0,0,147,,0,0,PGGFldElt,,0,0,PGGGloMod_Fld,,-25,-38,-38,-38,-38,-38
S,GlobalModel,A global model for F,0,1,0,0,0,0,0,0,0,PGGFldWrap,,PGGGloMod_Fld,-38,-38,-38,-38,-38
S,GlobalModel,A global model for F,0,1,0,0,0,0,0,0,0,PGGFldGrp,,PGGGloMod_Fld,-38,-38,-38,-38,-38
S,GlobalModel,"Global model with the given base model, and the corresponding hom W -> W2 where Gal(f)<W and Gal(model)<W2",0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global_Model,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,"Global model with the given base model, and the corresponding hom W -> W2 where Gal(f)<W and Gal(model)<W2",0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global_Model_Cheat,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,"Global model with the given base model, and the corresponding hom W -> W2 where Gal(f)<W and Gal(model)<W2",0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPolWrap,,0,0,PGGAlg_ResEval_Global_Model_Symmetric,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,"Global model with the given base model, and the corresponding hom W -> W2 where Gal(f)<W and Gal(model)<W2",0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPolGrp,,0,0,PGGAlg_ResEval_Global_Model_Symmetric,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,"Global model with the given base model, and the corresponding hom W -> W2 where Gal(f)<W and Gal(model)<W2",0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global_Model_Factors,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,Tower,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGFld,,0,0,PGGFld,,0,0,PGGAlg_ResEval_Global_Model_Tower,,82,-38,-38,-38,-38,-38
S,Tower,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGFld,,0,0,PGGFld,,0,0,PGGAlg_ResEval_Global_Model_RamTower,,82,-38,-38,-38,-38,-38
S,GlobalModel,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global_Model_Tower,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPol,,0,0,PGGAlg_ResEval_Global_Model_Select,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPolWrap,,0,0,PGGAlg_ResEval_Global_Model_RootOfUnity,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,GlobalModel,The tower from K up to L,0,3,0,0,0,0,0,0,0,PGGGloMod_Fld,,0,0,PGGPolWrap,,0,0,PGGAlg_ResEval_Global_Model_RootOfUniformizer,,PGGGloMod,PGGHomGrpPerm,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Rational,,82,-38,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Symmetric,,82,-38,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_RootOfUnity,,82,-38,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_RootOfUniformizer,,82,-38,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_PthRoots,,82,-38,-38,-38,-38,-38
S,ComplexEmbeddings,The complex embeddings of m extending embs,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Tower,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Symmetric,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_RootOfUnity,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_RootOfUniformizer,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_PthRoots,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Factors,,82,-38,-38,-38,-38,-38
S,ComplexRoots,The complex roots of the polynomial model,1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,PGGGloMod_Tower,,82,-38,-38,-38,-38,-38
S,AllComplexEmbeddings,A complex embedding into C,0,2,0,0,0,0,0,0,0,173,,0,0,PGGGloMod_Rational,,82,-38,-38,-38,-38,-38
S,AllComplexEmbeddings,A complex embedding into C,0,2,0,0,0,0,0,0,0,173,,0,0,PGGGloMod_Ext,,82,-38,-38,-38,-38,-38
