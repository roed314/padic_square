174,0
T,PGGArgSpec,-,0
A,PGGArgSpec,1,is_valid
T,PGGAttr_ArgSpecSubs,-,1,PGGAttr
A,PGGAttr_ArgSpecSubs,1,subs
S,PGG_Parsable,Tries to parse x according to A,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGArgSpec,,36,-1,-38,-38,-38,-38
S,PGG_Parsable,Tries to parse x according to A,0,2,0,0,0,0,0,0,0,PGGAttr_ArgSpecSubs,,0,0,PGGArgSpec,,36,-1,-38,-38,-38,-38
S,PGG_Parse,Parses x according to A,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGArgSpec,,-1,-38,-38,-38,-38,-38
S,PGG_ArgSpec,The argspec with the given is_valid function,0,1,0,0,0,0,0,0,0,41,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Union,The argspec for the union of the given choices,1,0,1,82,0,PGGArgSpec,1,0,0,0,0,0,0,0,82,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Transform,The argspec A with its output modified by tr,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGArgSpec,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_ISA,The argspec testing if the input has the given category,0,1,0,0,0,0,0,0,0,-1,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Literal,The argspec matching only X exactly,0,1,0,0,0,0,0,0,0,-1,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Attr,The argspec matching a PGGAttr,0,0,0,0,0,0,0,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Attr_NoArgs,The argspec matching a PGGAttr with no args. Outputs the name,0,0,0,0,0,0,0,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Attr,The argspec matching the PGGAttr with the given name and no args,0,1,0,0,0,0,0,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_AttrList,The argspec matching the PGGAttr with the given name and multiple args matching arg,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_AttrSelect,"The argspec matching the PGGAttr with the given name, and arguments [parg,rarg,parg,rarg,...,parg,rarg] (the final parg is optional). The pargs are predicates, and the rargs are results, so this is like select-else-select-else",0,4,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_AttrForEach,"The argspec matching a PGGAttr with the given name and three arguments [var,iter,arg] where var specifies a variable (or variables), iter is a list of items, each getting assigned to var and substituted into arg in turn. The output is the list of args",0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_AttrScope,"Parses `name[var,inner]` or `name[inner]` (with `var=dfltfar`), substrituting `var` for `mkval(var,inner)` in `inner`",0,4,0,0,0,0,0,0,0,PGGArgSpec,,0,0,41,,0,0,298,,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_AttrInt,The argspec matching the PGGAttr when it looks like an integer,0,0,0,0,0,0,0,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Attr,The argspec matching the PGGAttr with the given name and args,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,298,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Recursive,Allows the creation of a PGGArgSpec recursively,0,1,0,0,0,0,0,0,0,41,,PGGArgSpec,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Recursive,Allows the creation of n PGGArgSpecs recursively,0,2,0,0,0,0,0,0,0,41,,0,0,148,,82,-38,-38,-38,-38,-38
S,PGG_ArgSpec_ManyRecursive,Allows the creation of many PGGArgSpecs,0,1,0,0,0,0,0,0,0,41,,175,-38,-38,-38,-38,-38
S,PGG_ArgSpec_Predicate,The PGGArgSpec passing with the given test,0,1,0,0,0,0,0,0,0,41,,PGGArgSpec,-38,-38,-38,-38,-38
