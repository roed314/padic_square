174,0
T,PGGAlg_Stream,-,1,PGGAlg
T,PGGAlgState_Stream,-,1,PGGAlgState
A,PGGAlgState_Stream,1,overgroup
T,PGGAlg_Stream_Index,-,1,PGGAlg_Stream
A,PGGAlg_Stream_Index,3,filter,sort_key,dedupe
T,PGGAlgState_Stream_Index,-,1,PGGAlgState_Stream
A,PGGAlgState_Stream_Index,5,all_indices,indices,ii,subgroup_classes,stream_done
S,PGGAlg_Stream_Index_Make,"The ""Index"" streams algorithm",0,0,0,0,0,0,0,PGGAlg_Stream_Index,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_Stream_Index,,-38,-38,-38,-38,-38,-38
S,Start,Starts the algorithm,0,2,0,0,0,0,0,0,0,PGGGrpPerm,,0,0,PGGAlg_Stream_Index,,PGGAlgState_Stream_Index,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Stream,,0,0,PGGAlg_Stream,,PGGAlgState_Stream,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_Stream_Index,,0,0,PGGAlg_Stream_Index,,PGGAlgState_Stream_Index,-38,-38,-38,-38,-38
S,Reset,Resets s to its initial state,0,1,0,0,1,0,0,0,0,PGGAlgState_Stream_Index,,-38,-38,-38,-38,-38,-38
S,HasNextGroup,"True if the current stream has a next item. If so, also returns the item",0,1,0,0,0,0,0,0,0,PGGAlgState_Stream_Index,,36,-1,-38,-38,-38,-38
S,HasNextStream,"True if s has a next stream. If so, moves on to it",0,1,0,0,0,0,0,0,0,PGGAlgState_Stream_Index,,36,-38,-38,-38,-38,-38
