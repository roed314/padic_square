174,0
S,PGG_all_partition_groupings,All sequences Ps of length #Q of multisets of integers such that &join Ps eq P and &+Ps[i] eq Qs[i],2,0,1,198,0,148,1,1,82,0,148,2,0,0,0,0,0,0,0,82,,0,0,198,,82,-38,-38,-38,-38,-38
S,PGG_all_binnings,"If items is a sequence of multiplicities of different items, and bins is a sequence of multiplicities of different bins, returns a sequence of possible binnings, where a binning is a sequence (of distinct bins) of multisets (of bins) of multisets (of items) of valid binnings of the items. A multiset b of items into the ith bin is valid if is_valid(i,b) is true. If b is a partial binning extendable to a valid binning, then is_semivalid(i,b) must be true; this is used to terminate branches of the algorithm early",2,0,1,82,0,148,1,1,82,0,148,2,0,0,0,0,0,0,0,82,,0,0,82,,82,-38,-38,-38,-38,-38
S,PGG_has_random_subgroup,"Tries to find a subgroup H of G such that `is_valid(H)` is true. For all subgroups K of all such H, we must have `is_semivalid(K)` true. On success, returns true and a random such H. Otherwise, returns false. Can fail if MaxTries is exceeded",0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-27,,36,-27,-38,-38,-38,-38
S,PGG_has_random_subgroup_of_index,"Tries to find a subgroup H of G such that `is_valid(H)` is true. For all subgroups K of all such H, we must have `is_semivalid(K)` true. On success, returns true and a random such H. Otherwise, returns false. Can fail if MaxTries is exceeded",0,2,0,0,0,0,0,0,0,148,,0,0,-27,,36,-27,-38,-38,-38,-38
S,PGG_linear_divisions,"Given an interval of given length, find all possible divisions of it from among the given lengths. A division is a sequence of integers (sorted largest first) summing to the given length, whose elements are a subset of lengths. If Bound is given, the divisions must be lexicographically no more than it (i.e. the first element to disagree must be smaller)",1,1,1,198,0,148,2,0,0,0,0,0,0,0,198,,0,0,148,,82,-38,-38,-38,-38,-38
S,PGG_rectangle_divisions,"Given a rectangle of given width and height, returns all possible integer divisions of the rectangle with the given areas. A division is a series of vertical cuts, through the whole rectangle, followed by a series of horizontal cuts within each resulting rectangle. It is returned in the form `[<w1,[h11,h12,...]>,<w2,[h21,...]>,...]` where the `wi` sum to width, and for each `i` the `hij` sum to height, and areas is the multiset of `wi*hij`",1,2,1,198,0,148,3,0,0,0,0,0,0,0,198,,0,0,148,,0,0,148,,82,-38,-38,-38,-38,-38
S,PGG_realcomplexpairs,"Returns the sequence of `i` such that `xs[i]` are real, and the sequence of pairs `<i,j>` such that `xs[i]` and `xs[j]` are a complex conjugate pair. Assumes that xs are all the roots of a squarefree real polynomial in some order",1,0,1,82,0,172,1,0,0,0,0,0,0,0,82,,82,82,-38,-38,-38,-38
S,PGG_CloseVector,A vector of L close to vec,0,2,0,0,0,0,0,0,0,-1,,0,0,164,,165,-38,-38,-38,-38,-38