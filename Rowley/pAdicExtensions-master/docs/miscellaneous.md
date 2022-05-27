# Miscellaneous
{:#miscellaneous}

<a id="Evaluate"></a><a id="Evaluate--NwtnPgon--etc"></a><a id="Evaluate--NwtnPgon--FldRatElt"></a><a id="Evaluate--NwtnPgon--RngIntElt"></a>
> **Evaluate** (P :: *NwtnPgon*, x :: *FldRatElt*)
> 
> **Evaluate** (P :: *NwtnPgon*, x :: *RngIntElt*)
> 
> -> *FldRatElt*
> {:.ret}
{:.intrinsic}

Treating `P` as a piecewise-linear function, evaluates it at `x`.




<a id="LeftEvaluate"></a><a id="LeftEvaluate--NwtnPgon--etc"></a><a id="LeftEvaluate--NwtnPgon--FldRatElt"></a><a id="LeftEvaluate--NwtnPgon--RngIntElt"></a>
> **LeftEvaluate** (P :: *NwtnPgon*, x :: *FldRatElt*)
> 
> **LeftEvaluate** (P :: *NwtnPgon*, x :: *RngIntElt*)
> 
> -> *FldRatElt*
> {:.ret}
{:.intrinsic}

Evaluates the face of `P` to the left of `x` at `x`.




<a id="FactorialValuation"></a><a id="FactorialValuation--RngIntElt--etc"></a><a id="FactorialValuation--RngIntElt--RngIntElt"></a>
> **FactorialValuation** (n :: *RngIntElt*, p :: *RngIntElt*)
> 
> -> *RngIntElt*
> {:.ret}
{:.intrinsic}

The `p`-adic valuation of `n!` (equivalent to `Valuation(Factorial(n),p)` but more efficient).


<a id="BinomialValuation"></a><a id="BinomialValuation--RngIntElt--etc"></a><a id="BinomialValuation--RngIntElt--RngIntElt--RngIntElt"></a>
> **BinomialValuation** (n :: *RngIntElt*, k :: *RngIntElt*, p :: *RngIntElt*)
> 
> -> *RngIntElt*
> {:.ret}
{:.intrinsic}

The `p`-adic valuation of `n` choose `k`.


<a id="UnitFactorial"></a><a id="UnitFactorial--FldFin--etc"></a><a id="UnitFactorial--FldFin--RngIntElt"></a>
> **UnitFactorial** (F :: *FldFin*, n :: *RngIntElt*)
> 
> -> *FldFinElt*
> {:.ret}
{:.intrinsic}

The product of the integers up to `n` which are units in `F`. By Wilson's formula, if `n=kp+r` this is `(-1)^k * r!`.


<a id="ShiftedFactorial"></a><a id="ShiftedFactorial--FldFin--etc"></a><a id="ShiftedFactorial--FldFin--RngIntElt"></a>
> **ShiftedFactorial** (F :: *FldFin*, n :: *RngIntElt*)
> 
> -> *FldFinElt*
> {:.ret}
{:.intrinsic}

The product of the integers up to `n` shifted down to be units of `F`. This is the product of `UnitFactorial(n div p^i)` for all `i`.


<a id="ShiftedBinomial"></a><a id="ShiftedBinomial--FldFin--etc"></a><a id="ShiftedBinomial--FldFin--RngIntElt--RngIntElt"></a>
> **ShiftedBinomial** (F :: *FldFin*, n :: *RngIntElt*, k :: *RngIntElt*)
> 
> -> *FldFinElt*
> {:.ret}
{:.intrinsic}

The binomial coefficient `n` choose `k` shifted down to be a unit of `F`.


