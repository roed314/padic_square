174,0
T,PGGTimer,-,0
T,PGGTimerStack,-,1,PGGTimer
A,PGGTimerStack,1,current_label
T,PGGTimerBase,-,1,PGGTimerStack
A,PGGTimerBase,5,last_time,times,counts,labels,idxs
T,PGGTimerFork,-,1,PGGTimerStack
A,PGGTimerFork,2,parent,parent_label
T,PGGTimerNull,-,1,PGGTimer
V,PGG_Timer,2
S,PGG_Timer,A new timer,0,0,0,0,0,0,0,PGGTimerBase,-38,-38,-38,-38,-38
S,PGG_NullTimer,The null timer,0,0,0,0,0,0,0,PGGTimerNull,-38,-38,-38,-38,-38
S,Fork,Forks t,0,1,0,0,0,0,0,0,0,PGGTimerStack,,PGGTimerFork,-38,-38,-38,-38,-38
S,Fork,Forks t,0,1,0,0,0,0,0,0,0,PGGTimerNull,,PGGTimerNull,-38,-38,-38,-38,-38
S,FindIndex,"The index of the given label. If it doesn't exist, creates it first",1,1,1,82,0,298,2,0,0,0,0,0,0,0,82,,0,0,PGGTimerBase,,148,-38,-38,-38,-38,-38
S,Update,Assigns the accumulated time to the given label,1,2,1,82,0,298,3,0,0,1,0,0,0,0,82,,0,0,148,,0,0,PGGTimerBase,,-38,-38,-38,-38,-38,-38
S,Update,Assigns the accumulated time to the given label,1,2,1,82,0,298,3,0,0,1,0,0,0,0,82,,0,0,148,,0,0,PGGTimerFork,,-38,-38,-38,-38,-38,-38
S,Update,Assigns the accumulated time to the given label,1,2,1,82,0,298,3,0,0,1,0,0,0,0,82,,0,0,148,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Update,Assigns the accumulated time to the given label,1,1,1,82,0,298,2,0,0,1,0,0,0,0,82,,0,0,PGGTimer,,-38,-38,-38,-38,-38,-38
S,GetLabel,Gets the current label,0,1,0,0,0,0,0,0,0,PGGTimer,,82,-38,-38,-38,-38,-38
S,SetLabel,Sets the label,1,1,1,82,0,298,2,0,0,1,0,0,0,0,82,,0,0,PGGTimer,,-38,-38,-38,-38,-38,-38
S,Push,"Assings time to the current label, then appends x to the label",0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerBase,,-38,-38,-38,-38,-38,-38
S,Push,"Assings time to the current label, then appends x to the label",0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops n items from the label",0,2,0,0,1,0,0,0,0,148,,0,0,PGGTimerStack,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops n items from the label",0,2,0,0,1,0,0,0,0,148,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops n items from the label, the nth being x",0,3,0,0,1,0,0,0,0,298,,0,0,148,,0,0,PGGTimerStack,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops n items from the label, the nth being x",0,3,0,0,1,0,0,0,0,298,,0,0,148,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops the last item from the label, which must be x",0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimer,,-38,-38,-38,-38,-38,-38
S,Pop,"Assigns time to the current label, then pops the last item from the label",0,1,0,0,1,0,0,0,0,PGGTimer,,-38,-38,-38,-38,-38,-38
S,Swap,"Assigns time to the current label, then swaps the last item from the label to x",0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerStack,,-38,-38,-38,-38,-38,-38
S,Swap,"Assigns time to the current label, then swaps the last item from the label to x",0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Log,Assigns time to the current label with x appended,0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerStack,,-38,-38,-38,-38,-38,-38
S,Log,Assigns time to the current label with x appended,0,2,0,0,1,0,0,0,0,298,,0,0,PGGTimerNull,,-38,-38,-38,-38,-38,-38
S,Tree,A tree representation of t,0,1,0,0,0,0,0,0,0,PGGTimerBase,,303,-38,-38,-38,-38,-38
S,PrintTree,Prints t as a tree,0,1,0,0,1,0,0,0,0,PGGTimer,,-38,-38,-38,-38,-38,-38
