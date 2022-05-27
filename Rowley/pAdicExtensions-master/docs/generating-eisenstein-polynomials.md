# Generating Eisenstein polynomials
{:#generating-eisenstein-polynomials}



**Contents**
* [Creation of templates](#creation-of-templates)
* [Creation of polynomials](#creation-of-polynomials)

## Creation of templates
{:#creation-of-templates}

<a id="TemplateForEisensteinPolynomials"></a><a id="TemplateForEisensteinPolynomials--FldPadTmpl--etc"></a><a id="TemplateForEisensteinPolynomials--FldPadTmpl--RngIntElt"></a>
> **TemplateForEisensteinPolynomials** (F :: *FldPadTmpl*, n :: *RngIntElt*)
> 
> -> *Tmpl_RngUPolElt_FldPadTmpl*
> {:.ret}
{:.intrinsic}

Template for monic Eisenstein polynomials over `F` of degree `n`.

**Parameters**
- `J`: When given, restrict the discriminant valuation to `n+J-1`.

<a id="TemplateForEisensteinPolynomials-2"></a><a id="TemplateForEisensteinPolynomials--PadRamifPgon"></a>
> **TemplateForEisensteinPolynomials** (P :: *PadRamifPgon*)
> 
> -> *Tmpl_RngUPolElt_FldPadTmpl*
> {:.ret}
{:.intrinsic}

A template for Eisenstein polynomials with this ramification polygon.

**Parameters**
- `Limit := "Best"`: When "Best", make the template as small as possible; when "Disc" use the discriminant bound; otherwise it must be an integer.
- `Residues := true`: When true, incorporate residue information when it is available

## Creation of polynomials
{:#creation-of-polynomials}

<a id="#"></a><a id="#--Tmpl_FldPadTmplElt"></a><a id="#--Tmpl_RngUPolElt_FldPadTmpl"></a>
> **\'#\'** (T :: *Tmpl_FldPadTmplElt*)
> 
> **\'#\'** (T :: *Tmpl_RngUPolElt_FldPadTmpl*)
> 
> -> *RngIntElt*
> {:.ret}
{:.intrinsic}

The number of elements in the template.




<a id="ToSequence"></a><a id="ToSequence--Tmpl_FldPadTmplElt"></a><a id="ToSequence--Tmpl_RngUPolElt_FldPadTmpl"></a>
> **ToSequence** (T :: *Tmpl_FldPadTmplElt*)
> 
> **ToSequence** (T :: *Tmpl_RngUPolElt_FldPadTmpl*)
> 
> -> []
> {:.ret}
{:.intrinsic}

The elements of `T` as a sequence.




<a id="Random"></a><a id="Random--Tmpl_FldPadTmplElt"></a><a id="Random--Tmpl_RngUPolElt_FldPadTmpl"></a>
> **Random** (T :: *Tmpl_FldPadTmplElt*)
> 
> -> *FldPadTmplElt*
> {:.ret}
> 
> **Random** (T :: *Tmpl_RngUPolElt_FldPadTmpl*)
> 
> -> *RngUPolElt_FldPadTmpl*
> {:.ret}
{:.intrinsic}

A random element of `T`.




<a id="IndexSet"></a><a id="IndexSet--Tmpl_FldPadTmplElt"></a><a id="IndexSet--Tmpl_RngUPolElt_FldPadTmpl"></a>
> **IndexSet** (T :: *Tmpl_FldPadTmplElt*)
> 
> **IndexSet** (T :: *Tmpl_RngUPolElt_FldPadTmpl*)
> 
> -> Any
> {:.ret}
{:.intrinsic}

A set of indices for elements of `T`. This is quick to generate and can be iterated over.




<a id="GetIndex"></a><a id="GetIndex--Tmpl_FldPadTmplElt--etc"></a><a id="GetIndex--Tmpl_FldPadTmplElt--any"></a><a id="GetIndex--Tmpl_RngUPolElt_FldPadTmpl--etc"></a><a id="GetIndex--Tmpl_RngUPolElt_FldPadTmpl--any"></a>
> **GetIndex** (T :: *Tmpl_FldPadTmplElt*, idx)
> 
> -> *FldPadTmplElt*
> {:.ret}
> 
> **GetIndex** (T :: *Tmpl_RngUPolElt_FldPadTmpl*, idx)
> 
> -> *RngUPolElt_FldPadTmpl*
> {:.ret}
{:.intrinsic}

Gets the element of `T` at index `idx` (which should be an element of the index set).




