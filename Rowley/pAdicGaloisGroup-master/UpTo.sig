174,0
T,PGGAlg_GaloisGroup_ARM_UpTo,-,1,PGGAlg
T,PGGAlg_GaloisGroup_ARM_UpTo_HomConj,-,1,PGGAlg_GaloisGroup_ARM_UpTo
T,PGGAlg_GaloisGroup_ARM_UpTo_Full,-,1,PGGAlg_GaloisGroup_ARM_UpTo_HomConj
T,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding,-,1,PGGAlg_GaloisGroup_ARM_UpTo_HomConj
A,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding,1,check_injective
T,PGGAlg_GaloisGroup_ARM_UpTo_Symmetric,-,1,PGGAlg_GaloisGroup_ARM_UpTo_HomConj
T,PGGAlgState_GaloisGroup_ARM_UpTo,-,1,PGGAlgState
T,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj,-,1,PGGAlgState_GaloisGroup_ARM_UpTo
A,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj,2,hom,group
T,PGGAlgState_GaloisGroup_ARM_UpTo_Symmetric,-,1,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj
T,PGGAlgState_GaloisGroup_ARM_UpTo_Full,-,1,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj
T,PGGAlgState_GaloisGroup_ARM_UpTo_ResolventEmbedding,-,1,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj
S,PGGAlg_GaloisGroup_ARM_UpTo_Symmetric_Make,Used to signify that groups are determined up to the symmetric group only,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Symmetric,-38,-38,-38,-38,-38
S,PGGAlg_GaloisGroup_ARM_UpTo_Full_Make,Used to signify that the Galois group is determined up to W-conjugacy,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Full,-38,-38,-38,-38,-38
S,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding_Make,Used to signify that the Galois group is determined up to Wtil-conjugacy under e:W->Wtil,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Full,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Symmetric,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding,,-38,-38,-38,-38,-38,-38
S,Start,Starts the algorithm,0,2,0,0,0,0,0,0,0,PGGHomGrpPerm,,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,PGGAlgState_GaloisGroup_ARM_UpTo,-38,-38,-38,-38,-38
S,Start,Starts the algorithm,0,2,0,0,0,0,0,0,0,PGGHomGrpPerm,,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Full,,PGGAlgState_GaloisGroup_ARM_UpTo_Full,-38,-38,-38,-38,-38
S,Start,Starts the algorithm,0,2,0,0,0,0,0,0,0,PGGHomGrpPerm,,0,0,PGGAlg_GaloisGroup_ARM_UpTo_Symmetric,,PGGAlgState_GaloisGroup_ARM_UpTo_Symmetric,-38,-38,-38,-38,-38
S,Start,Starts the algorithm,0,2,0,0,0,0,0,0,0,PGGHomGrpPerm,,0,0,PGGAlg_GaloisGroup_ARM_UpTo_ResolventEmbedding,,PGGAlgState_GaloisGroup_ARM_UpTo_ResolventEmbedding,-38,-38,-38,-38,-38
S,IsEquivalent,True if G1 and G2 are equivalent,0,3,0,0,0,0,0,0,0,224,,0,0,224,,0,0,PGGAlgState_GaloisGroup_ARM_UpTo,,36,-38,-38,-38,-38,-38
S,IsEquivalent,True if G1 and G2 are equivalent,0,3,0,0,0,0,0,0,0,224,,0,0,224,,0,0,PGGAlgState_GaloisGroup_ARM_UpTo_HomConj,,36,-38,-38,-38,-38,-38
