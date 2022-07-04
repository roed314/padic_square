Attach("getslopes.m");
Attach("padic_square/Montes/+IdealsNF.m");
load "polys23";
"Polynomials are in polys23";

AddAttribute(FldRat, "slopecache");
AddAttribute(FldRat, "galdat");
AddAttribute(FldRat, "dat16");

low:=StringToInteger(low);

for j:=low to low+19 do get_slopes(polys23[j],2); end for;
