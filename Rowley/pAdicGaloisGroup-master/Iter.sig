174,0
T,PGGIter,-,0
A,PGGIter,1,cache
T,PGGIter_List,-,1,PGGIter
A,PGGIter_List,2,list,idx
T,PGGIter_Filter,-,1,PGGIter
A,PGGIter_Filter,2,iter,filter
T,PGGIter_Truncate,-,1,PGGIter
A,PGGIter_Truncate,2,iter,limit
T,PGGIter_Enumerate,-,1,PGGIter
A,PGGIter_Enumerate,2,iter,count
T,PGGIter_Apply,-,1,PGGIter
A,PGGIter_Apply,2,iter,map
T,PGGIter_User,-,1,PGGIter
A,PGGIter_User,2,state,has_next
S,PGG_Iter,Makes an iterator X with whose state attribute is State and `HasNext(X)` is `has_next(X)`,0,1,0,0,0,0,0,0,0,41,,PGGIter,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter_User,,36,-1,-38,-38,-38,-38
S,PGG_ToIter,Makes an iterator,0,1,0,0,0,0,0,0,0,82,,PGGIter,-38,-38,-38,-38,-38
S,PGG_ToIter,Makes an iterator,0,1,0,0,0,0,0,0,0,168,,PGGIter,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter_List,,36,-1,-38,-38,-38,-38
S,ToList,Converts X to a list,0,1,0,0,0,0,0,0,0,PGGIter,,168,-38,-38,-38,-38,-38
S,ToList,Converts X to a list,0,1,0,0,0,0,0,0,0,PGGIter_List,,168,-38,-38,-38,-38,-38
S,ToSequence,Converts X to a sequence,0,1,0,0,0,0,0,0,0,PGGIter,,82,-38,-38,-38,-38,-38
S,Enumerate,"The sequence <i,x_i> for x_i in X",0,1,0,0,0,0,0,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter_Enumerate,,36,-1,-38,-38,-38,-38
S,Apply,The sequence f(x) for x in X,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter,,36,-1,-38,-38,-38,-38
S,Filter,The x in X such that f(x) is true,0,2,0,0,0,0,0,0,0,-1,,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter_Filter,,36,-1,-38,-38,-38,-38
S,Reverse,The reverse,0,1,0,0,0,0,0,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,Shuffle,Randomize the order of X,0,1,0,0,0,0,0,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,SortBy,Sorts X so that keyfunc(x) are naturally sorted,0,2,0,0,0,0,0,0,0,41,,0,0,PGGIter,,PGGIter,-38,-38,-38,-38,-38
S,FilterHasNext,"True if there is `x` in X such that `f(x)` is true. If so, returns the first such `x` and its index",0,2,0,0,0,0,0,0,0,-1,,0,0,PGGIter,,36,-1,148,-38,-38,-38
S,Truncate,Limits the length of X to n,0,2,0,0,0,0,0,0,0,148,,0,0,PGGIter,,PGGIter_Truncate,-38,-38,-38,-38,-38
S,Truncate,Limits the length of X to n,0,2,0,0,0,0,0,0,0,148,,0,0,PGGIter_List,,PGGIter_List,-38,-38,-38,-38,-38
S,HasNext,True if X has a next item,0,1,0,0,0,0,0,0,0,PGGIter,,36,-1,-38,-38,-38,-38
