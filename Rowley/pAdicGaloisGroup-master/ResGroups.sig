174,0
T,PGGAlg_ResGroups,-,1,PGGAlg
T,PGGAlgState_ResGroups,-,1,PGGAlgState
A,PGGAlgState_ResGroups,2,overgroup_embedding,upto
T,PGGAlg_ResGroups_All,-,1,PGGAlg_ResGroups
A,PGGAlg_ResGroups_All,4,subgroup_choice,statistic,is_subgroup_of,dedupe
T,PGGAlgState_ResGroups_All,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_All,4,galois_conjugacy,resolvent_overgroup,possible_groups,choice_state
T,PGGAlg_ResGroups_RootsMaximal,-,1,PGGAlg_ResGroups
A,PGGAlg_ResGroups_RootsMaximal,1,dedupe
T,PGGAlgState_ResGroups_RootsMaximal,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_RootsMaximal,5,galois_conjugacy,resolvent_overgroup,cur_group,possible_subgroups,conjugacy_classes
T,PGGAlg_ResGroups_Sequence,-,1,PGGAlg_ResGroups
A,PGGAlg_ResGroups_Sequence,1,seq
T,PGGAlgState_ResGroups_Sequence,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_Sequence,3,idx,cur_state,prev_state
T,PGGAlg_ResGroups_Null,-,1,PGGAlg_ResGroups
T,PGGAlgState_ResGroups_Null,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_Null,1,previous
T,PGGAlg_ResGroups_Maximal,-,1,PGGAlg_ResGroups
A,PGGAlg_ResGroups_Maximal,8,subgroup_choice,statistic,descend_when,descend,useful,reprocess,reset,dedupe
T,PGGAlgState_ResGroups_Maximal,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_Maximal,8,galois_conjugacy,resolvent_overgroup,choice_state,graph,pool,resolvents,useful_is_special,actually_useful_is_special
T,PGG_ResGroups_Graph,PGG_ResGroups_GraphNode,0
A,PGG_ResGroups_Graph,5,galois_conjugacy,dedupe,conjugacy_classes,id_to_node,class_to_id
A,PGG_ResGroups_GraphNode,9,parent,class,id,children,possible_children,parents,possible_parents,maybe_equal,maybe_subgroup
T,PGGAlgState_ResGroups_Maximal_Graph,PGGAlgState_ResGroups_Maximal_GraphNode,1,PGG_ResGroups_Graph
A,PGGAlgState_ResGroups_Maximal_Graph,1,state
T,PGGAlgState_ResGroups_Maximal_GraphNode,-,1,PGG_ResGroups_GraphNode
A,PGGAlgState_ResGroups_Maximal_GraphNode,1,pooled_children
T,PGGAlg_ResGroups_Maximal2,-,1,PGGAlg_ResGroups
A,PGGAlg_ResGroups_Maximal2,5,subgroup_choice,statistic,reset,dedupe,descend
T,PGGAlgState_ResGroups_Maximal2,-,1,PGGAlgState_ResGroups
A,PGGAlgState_ResGroups_Maximal2,5,galois_conjugacy,resolvent_overgroup,choice_state,graph,pool
S,PGGAlg_ResGroups_All_Make,"The ""All"" resolvent groups algorithm",0,0,0,0,0,0,0,PGGAlg_ResGroups_All,-38,-38,-38,-38,-38
S,PGGAlg_ResGroups_RootsMaximal_Make,"The ""RootsMaximal"" resolvent groups algorithm",0,0,0,0,0,0,0,PGGAlg_ResGroups_RootsMaximal,-38,-38,-38,-38,-38
S,PGGAlg_ResGroups_Sequence_Make,"The ""Sequence"" resolvent groups algorithm",0,1,0,0,0,0,0,0,0,-1,,PGGAlg_ResGroups_Sequence,-38,-38,-38,-38,-38
S,PGGAlg_ResGroups_Null_Make,"The ""Null"" resolvent groups algorithm",0,0,0,0,0,0,0,PGGAlg_ResGroups_Null,-38,-38,-38,-38,-38
S,PGGAlg_ResGroups_Maximal_Make,"The ""Maximal"" resolvent groups algorithm",0,0,0,0,0,0,0,PGGAlg_ResGroups_Maximal,-38,-38,-38,-38,-38
S,PGGAlg_ResGroups_Maximal2_Make,"The ""Maximal2"" resolvent groups algorithm",0,0,0,0,0,0,0,PGGAlg_ResGroups_Maximal2,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_All,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_RootsMaximal,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_Sequence,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_Null,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_Maximal,,-38,-38,-38,-38,-38,-38
S,Print,Print,0,1,0,0,1,0,0,0,0,PGGAlg_ResGroups_Maximal2,,-38,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_All,,PGGAlgState_ResGroups_All,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_RootsMaximal,,PGGAlgState_ResGroups_RootsMaximal,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_Sequence,,PGGAlgState_ResGroups_Sequence,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_Null,,PGGAlgState_ResGroups_Null,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_Maximal,,PGGAlgState_ResGroups_Maximal,-38,-38,-38,-38,-38
S,Start,Starts the algorithm and returns its state,0,3,0,0,0,0,0,0,0,PGGAlg_GaloisGroup_ARM_UpTo,,0,0,PGGHomGrpPerm,,0,0,PGGAlg_ResGroups_Maximal2,,PGGAlgState_ResGroups_Maximal2,-38,-38,-38,-38,-38
S,RestartFromScratch,"Same as ``Start(alg, s`overgroup_embedding)`` but with a warning that it's happening",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResGroups,,PGGAlgState_ResGroups,-38,-38,-38,-38,-38
S,AllowRestartFromScratch,True if we allow restarting an instance of b from a,0,2,0,0,0,0,0,0,0,PGGAlg_ResGroups,,0,0,PGGAlg_ResGroups,,36,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResGroups,,PGGAlgState_ResGroups,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlg_ResGroups_Null,,PGGAlgState_ResGroups_Null,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Null,,0,0,PGGAlg_ResGroups_Null,,PGGAlgState_ResGroups_Null,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Null,,0,0,PGGAlg_ResGroups,,PGGAlgState_ResGroups,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_All,,0,0,PGGAlg_ResGroups_All,,PGGAlgState_ResGroups_All,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal2,,0,0,PGGAlg_ResGroups_Maximal2,,PGGAlgState_ResGroups_Maximal2,-38,-38,-38,-38,-38
S,Restart,"Like ``Start(alg, s`overgroup_embedding)`` but tries to preserve information",0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,0,0,PGGAlg_ResGroups_Sequence,,PGGAlgState_ResGroups_Sequence,-38,-38,-38,-38,-38
S,Init,Initializes X,0,3,0,0,1,0,0,0,0,PGGAlg_SubgrpDedupe,,0,0,PGGConj,,0,0,PGG_ResGroups_Graph,,-38,-38,-38,-38,-38,-38
S,Init,Initializes X,0,2,0,0,1,0,0,0,0,PGGAlgState_ResGroups,,0,0,PGGAlgState_ResGroups_Maximal_Graph,,-38,-38,-38,-38,-38,-38
S,PGGAlgState_ResGroups_Maximal_Graph_Make,The graph for s,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,PGGAlgState_ResGroups_Maximal_Graph,-38,-38,-38,-38,-38
S,PGGAlgState_ResGroups_Maximal_Graph_Make,The graph for s,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal2,,PGGAlgState_ResGroups_Maximal_Graph,-38,-38,-38,-38,-38
S,Init,Initializes n,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,Init,Initializes n,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-38,-38,-38,-38,-38,-38
S,Node,The node of X at G,0,2,0,0,0,0,0,0,0,_PGGSubgrpcls,,0,0,PGG_ResGroups_Graph,,PGG_ResGroups_GraphNode,-38,-38,-38,-38,-38
S,Node,The node of X at G,0,2,0,0,0,0,0,0,0,224,,0,0,PGG_ResGroups_Graph,,PGG_ResGroups_GraphNode,-38,-38,-38,-38,-38
S,Node,The node of X at G,0,2,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,0,0,PGG_ResGroups_Graph,,PGG_ResGroups_GraphNode,-38,-38,-38,-38,-38
S,IsCoercible,True if G is coercible to a node of X,0,2,0,0,0,0,0,0,0,-1,,0,0,PGG_ResGroups_Graph,,36,PGG_ResGroups_Graph,-38,-38,-38,-38
S,Root,The root node of X,0,1,0,0,0,0,0,0,0,PGG_ResGroups_Graph,,PGG_ResGroups_GraphNode,-38,-38,-38,-38,-38
S,Parent,Parent,0,1,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,PGG_ResGroups_Graph,-38,-38,-38,-38,-38
S,Hash,Hash,0,1,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,148,-38,-38,-38,-38,-38
S,eq,Equality,0,2,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,0,0,PGG_ResGroups_GraphNode,,36,-38,-38,-38,-38,-38
S,Children,The child nodes of n,0,1,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,83,-38,-38,-38,-38,-38
S,PossibleChildren,The possible children of n,0,1,0,0,0,0,0,0,0,PGG_ResGroups_GraphNode,,83,-38,-38,-38,-38,-38
S,PossibleChildrenHook,Overload this to do something when `PossibleChildren(n)` is called for the first time,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,PossibleChildrenHook,Overload this to do something when `PossibleChildren(n)` is called for the first time,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotEqual,Record that n is not equal to the Galois group,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotEqualHook,Overload this to do something when `NotEqual(n)` is called,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotEqualHook,Overload this to do something when `NotEqual(n)` is called,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotSubgroup,Record that the Galois group is not a subgroup of n,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotSubgroupHook,Overload this to do something when `NotSubgroup(n)` is called,0,1,0,0,1,0,0,0,0,PGG_ResGroups_GraphNode,,-38,-38,-38,-38,-38,-38
S,NotSubgroupHook,Overload this to do something when `NotSubgroup(n)` is called,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-38,-38,-38,-38,-38,-38
S,PGG_Image,The image of G under h,0,2,0,0,0,0,0,0,0,175,,0,0,224,,224,-38,-38,-38,-38,-38
S,PGG_Image,"The image of G under hs, the corresponding GSets, and the combined map",1,1,1,82,0,175,2,0,0,0,0,0,0,0,82,,0,0,224,,224,82,175,-38,-38,-38
S,GroupStat,"The image of G under hs, the corresponding GSets, and the combined map",0,3,0,0,0,0,0,0,0,175,,0,0,224,,0,0,PGGStat,,-1,-38,-38,-38,-38,-38
S,GroupStat,"The image of G under hs, the corresponding GSets, and the combined map",1,2,1,82,0,175,3,0,0,0,0,0,0,0,82,,0,0,224,,0,0,PGGStat,,-1,-38,-38,-38,-38,-38
S,GroupStat,The statistic of n under the action h,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-1,-38,-38,-38,-38,-38
S,CurrentState,The current state,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,PGGAlgState_ResGroups,-38,-38,-38,-38,-38
S,NextState,"The next state. At the end of the sequence, returns Null",0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,PGGAlgState_ResGroups,-38,-38,-38,-38,-38
S,AtEnd,True if s has reached the end of the sequence,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,36,-38,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_All,,36,224,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_RootsMaximal,,36,224,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,36,224,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Null,,36,224,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,36,224,-38,-38,-38,-38
S,IsDone,True if we have deduced the Galois group,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal2,,36,224,-38,-38,-38,-38
S,TheGroup,"The Galois group, assuming it is already known",0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,224,-38,-38,-38,-38,-38
S,PossibleSubgroups,The subgroups of s`cur_group which have not yet been ruled out,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_RootsMaximal,,82,-38,-38,-38,-38,-38
S,ProcessResolvent,Use the resolvent to get information about the Galois group,0,3,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,PGGAlgState_ResGroups_All,,-38,-38,-38,-38,-38,-38
S,ProcessResolvent,Use the resolvent to get information about the Galois group,0,3,0,0,1,0,0,0,0,224,,0,0,PGGPol,,0,0,PGGAlgState_ResGroups_RootsMaximal,,-38,-38,-38,-38,-38,-38
S,ProcessResolvent,Use the resolvent to get information about the Galois group,0,3,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,PGGAlgState_ResGroups_Sequence,,-38,-38,-38,-38,-38,-38
S,ProcessResolvent,Use the resolvent to get information about the Galois group,0,3,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,PGGAlgState_ResGroups_Maximal,,-38,-38,-38,-38,-38,-38
S,Stabilizers,The preimages under h of the stabilizers of v in G@h,0,3,0,0,0,0,0,0,0,PGGStatVal,,0,0,175,,0,0,224,,82,-38,-38,-38,-38,-38
S,Stabilizers,The preimages under h of the stabilizers of v in G@h,1,1,1,82,0,175,3,0,0,0,0,0,0,0,PGGStatVal,,0,0,82,,0,0,224,,82,-38,-38,-38,-38,-38
S,ProcessResolvent,The preimages under h of the stabilizers of v in G@h,0,3,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,PGGAlgState_ResGroups_Maximal2,,-38,-38,-38,-38,-38,-38
S,ProcessResolvent,Processes the ridxth resolvent of s,0,2,0,0,1,0,0,0,0,148,,0,0,PGGAlgState_ResGroups_Maximal,,-38,-38,-38,-38,-38,-38
S,MustDescend,True if we must descend,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,36,-38,-38,-38,-38,-38
S,CanDescend,True if we can descend,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,36,-38,-38,-38,-38,-38
S,ShouldDescend,True if we should descend,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,36,-38,-38,-38,-38,-38
S,PooledChildren,The pooled children of n,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,83,-38,-38,-38,-38,-38
S,PoolableChildren,The poolable children of n,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,83,-38,-38,-38,-38,-38
S,RemoveFromPool,Removes n from the pool,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal_GraphNode,,-38,-38,-38,-38,-38,-38
S,AddToPool,Adds the given nodes to the pool,1,0,1,83,0,PGGAlgState_ResGroups_Maximal_GraphNode,1,0,0,1,0,0,0,0,83,,-38,-38,-38,-38,-38,-38
S,Descend,Descends,0,1,0,0,1,0,0,0,0,PGGAlgState_ResGroups_Maximal,,-38,-38,-38,-38,-38,-38
S,UsefulHook,This is called when i is useful. Overload this,0,2,0,0,1,0,0,0,0,PGGAlgState_TrancheItem,,0,0,PGGAlgState_ResGroups,,-38,-38,-38,-38,-38,-38
S,UsefulHook,This is called when i is useful. Overload this,0,2,0,0,1,0,0,0,0,PGGAlgState_TrancheItem,,0,0,PGGAlgState_ResGroups_Maximal,,-38,-38,-38,-38,-38,-38
S,IsUseful,True if i represents a useful group,0,2,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,0,0,PGGAlgState_ResGroups,,36,-38,-38,-38,-38,-38
S,IsUseful,True if i represents a useful group,0,2,0,0,0,0,0,0,0,PGGAlgState_TrancheItem,,0,0,PGGAlgState_ResGroups_Maximal,,36,-38,-38,-38,-38,-38
S,IsUseful,True if G is useful,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGAlgState_ResGroups_All,,36,-38,-38,-38,-38,-38
S,IsUseful,True if G is useful,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGAlgState_ResGroups_Maximal,,36,-38,-38,-38,-38,-38
S,IsUseful,True if G is useful,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGAlgState_ResGroups_Maximal2,,36,-38,-38,-38,-38,-38
S,Information,The information provided by this item,0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_All,,0,0,PGGAlgState_TrancheItem,,402,-38,-38,-38,-38,-38
S,Diversity,The information provided by this item,0,2,0,0,0,0,0,0,0,PGGAlgState_ResGroups_All,,0,0,PGGAlgState_TrancheItem,,148,-38,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_All,,36,224,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_RootsMaximal,,36,224,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Sequence,,36,224,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Null,,36,224,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal,,36,224,-38,-38,-38,-38
S,HasSubgroup,True if s has a subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups_Maximal2,,36,224,-38,-38,-38,-38
S,Subgroup,A subgroup to form a resolvent from,0,1,0,0,0,0,0,0,0,PGGAlgState_ResGroups,,224,-38,-38,-38,-38,-38
