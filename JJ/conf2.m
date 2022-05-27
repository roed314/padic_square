Attach("getslopes.m");
load "polys23";
"Polynomials are in polys23";

low:=StringToInteger(low);

for j:=low to low+10 do get_slopes(polys23[j],2); end for;
