# To be done before imminent pull request

 * Update treatment of Swan vs Artin slopes to make it more uniform across field and family homepages and search pages by making the following changes:
   * In "Fields" section of family homepage, add "Visible Swan slopes" and "Swan slope content" columns, hidden by default.
 * Work on knowls
   * [lf.eisenstein_diagram](https://olive.lmfdb.xyz/knowledge/edit/lf.eisenstein_diagram) - Explaining all the features of the Eisenstein diagram and how it connects to the invariants of the family
   * [lf.family_ambiguity](https://olive.lmfdb.xyz/knowledge/edit/lf.family_ambiguity) and [lf.family_mass](https://olive.lmfdb.xyz/knowledge/edit/lf.family_mass) need updating (the first should be the definition in 3.5 for the ambiguity of I/K, and the second should include formulas for relative and absolute mass.
   * In [lf.slope_content](https://olive.lmfdb.xyz/knowledge/edit/lf.slope_content) we need to clarify whether we are talking about Swan or Artin slopes.
   * We need to update [lf.field.label](https://olive.lmfdb.xyz/knowledge/edit/lf.field.label) to describe the new labels, and labels for families and subfamilies.
   * Update [lf.visible_slopes](https://olive.lmfdb.xyz/knowledge/show/lf.visible_slopes) to link to [lf.hidden_slopes](https://olive.lmfdb.xyz/knowledge/show/lf.hidden_slopes).
   * Update completeness, reliability, and label knowls.
   * Clarify difference (in knowls) between mass and absolute mass in cases where f is not 1 (e.g. [here](https://olive.lmfdb.xyz/padicField/family/2.2.8.54a)).  Also relevant for mass vs mass stored columns in family search results.
 * Update other parts of the LMFDB (notably the number field pages) to use new labels -- if we transition the labels column to hold the new labels (copying the current label column to old_labels), then this will happen automatically
 * Update other parts of the LMFDB (notably the number field pages) to use new labels.  Make sure https://olive.lmfdb.xyz/NumberField/?completions=5.4.1.0a1.1 works (this is the link to "number fields with this as the completion) (JJ)
 * Switch ramification polygon style to Santi's, deal with points above the polygon (e.g. in 2.1.16.60f1.228) (DLR)
 * Deal with very close lines in the Eisenstein diagram (e.g. in 2.1.32.130bk) (DLR)
 * Make it possible to search for fields containing a given one
 * Refresh stats once everything computed, make sure tables go up to 16 (DLR)
 * Review knowls (DLR)
 * Tiny fix: On family page, when jump set is empty use "[\ ]" instead of "[]"

# Issues for later

## Front-end changes (most of these should be collected into Github issues)
 * Add the canonical subfields to search columns (on family page and field search results)
 * Check that the new hidden column is correctly sorted, update hidden_slopes attribute in family.py to use it.
 * Improve the disply in the Varying section of a family homepage in cases where there are many different Galois groups and hidden slopes (see [here](https://olive.lmfdb.xyz/padicField/family/2.2.8.54a) for an example).  This could take the form of making it a table of some sort, or adding information on degree, which is currently less apparent for Galois groups.
 * Search for families (and maybe fields) based on: slope multiplicity, slopes, heights, rams.
 * In the field search results on a family homepage, have columns for ai, bi, ci (better alignment, smaller)
 * After reloading data, fix mass and mass_stored display to use mixed fractions
 * Need ctr0 for base sorting on family, need sort counter for sorting by slopes (grouped by base and n?)
 * Add sort options for family searches (including sort keys on base for relative families)
 * Run PolRedPadic on the Find box so that users can type in any polynomial and get the corresponding p-adic field (echoing functionality for number fields)
 * Name tame fields as Q_q(nth_root(pi)) for an appropriate pi, then update the base method in family.py.  Should also name two important families of wildly ramified extensions: pth power roots of unity (subsumed by general cyclotomic fields), and pth power roots of any element in the tame subfield.
 * It would be nice to show actual automorphism group rather than just its order
 * Add ability to compute families on the fly (in larger degree for example)
 * Add more dynamic columns to dynamic stats, fix links (currently if you click on the 32 entry [here](https://olive.lmfdb.xyz/padicField/dynamic_stats?p=2&n=8&visible_quantifier=exactly&visible=[2,+3,+17%2F4]&col1=slopes&totals1=yes&col2=galois_label&totals2=yes&proportions=none) it also includes fields from another family with c=26 since visible isn't being unparsed correctly).
 * Add a toggle to display defining polynomial of a field in Eisenstein form vs just as powers of x (analogous to homogenize/simplify for elliptic curves)
 * Commands for Sage/Magma
 * On a field homepage, add the other canonical subfields and links to field pages for fields in the tower.


## LMFDB data

 * Compute the field labels for each polynomial in families and store in another table.
 * add relative defining polynomials and relative Galois groups to the database (Galois group over the canonical subfields).

## Other code

 * In family_generation.py, fix `__iter__` when base is not Qp.
 * Given a non-Galois K/Qp, write code to find the slope filtration on Gal(K/Qp) as a sequence of subgroups.  John suggested first finding the slopes (giving the sizes for the groups in the filtration), then trying to use resolvents to determine which options have the correct sized fixed field.
 * Given a non-Galois K/Qp, write code to find the canonical filtration by visible slopes Q_p < K_1 < K_2 < ... < K_j < K, where each extension K_j/K_{j-1} has a single slope associated to it, strictly increasing.  Hopefully this can run in moderate degree (e.g. 1000 over Q2), and may be easier than the general find-subfields problem.

## Theory pursuits

 * Investigate how Galois groups vary in a family (is there a generic group?  can we predict which groups will arise a priori somehow?)
  * Monge-like formula for number of fields in a family

## Paper

 * Write second paper

## Other

 * Review our previous reports and pick up threads that we haven't been working on.
 * David Roberts mentioned work in number fields for finding Galois fields with root discriminant bounded by a specific value (around 45).  This is connected to the compositum game: trying to find different number fields with very similar ramification so that their compositum has small root discriminant.  For this, it would be helpful to be able to search based on the following partial order on number fields: For every prime we associate two numbers: the top wild slope s (which will be 1 for tame and 0 for unramified) and the lcm t of the tame degrees above that prime.  We say K <= L if s(K) <= s(L) and s(K) | s(L) for all primes p.  Given a number field L we want to find all other K in the database that are less than or equal to L.

# Done

 * ~~Fix bug where variable numbering is repeated (e.g. [2-4.4_10_11_11](https://olive.lmfdb.xyz/padicField/family/2-4.4_10_11_11)).  Results from bands that have the same top, so the same red point.~~
 * ~~Fix bug where fractional heights lead to fractional variable subscripts (e.g. [2-6.6_12_13_13](https://olive.lmfdb.xyz/padicField/family/2-6.6_12_13_13))~~
 * ~~Fix bug where red and green points can overlap (e.g. [2-2.2_5_7_9](https://olive.lmfdb.xyz/padicField/family/2-2.2_5_7_9))~~
 * ~~Fix broken [random family](https://olive.lmfdb.xyz/padicField/families/?n=8&search_type=Random)~~
 * ~~Allow families to include a visible tame degree~~
 * ~~Figure out how to vary the base field~~
 * ~~When the denominator of a slope is not a power of p, the formula for green points, (u_i′,v_i) = (⟨h_i′⟩,⌈h_i′⌉) given on page 17 of David's notes breaks since the u-coordinate is not integral.  The number of bands has also decreased, so we no longer have the same number of points.  What's the right analogue of green points in this setting?  Do they just not exist, since there isn't an integral point at the bottom of the band (yes, we think that's right)?~~
 * ~~Add the size of the automorphism group to the rows in the list of fields at the bottom of each family page.~~
 * ~~Vertical scaling on family picture can be very bad (e.g. https://olive.lmfdb.xyz/padicField/family/3-2.1_1_1)~~
 * ~~Fix tight boundaries around picture, which are cutting off the top of a red diamond in /padicField/family/3-2.1_2~~
 * ~~Search by and display order of automorphism group~~
 * ~~Group by finer invariants in list of fields, with link to expand (like in ECQ isogeny classes)~~
 * ~~Fix mass~~
 * ~~Create lf_families table and switch families search to use it.~~
 * ~~Create a lf_families table to serve as a table behind the [families searches](https://olive.lmfdb.xyz/padicField/families/).  Make sure that we're happy with our labeling scheme first, and figure out how we want to handle unramified and tame bases for families.~~
 * ~~Run Keating/JumpSetPack.m at scale to compute jump sets across lf_fields.~~
 * ~~Think about how to reduce the number of variables in cases like [2-2.2_5_7_9](https://olive.lmfdb.xyz/padicField/family/2-2.2_5_7_9), or divide into subfamilies.~~
 * ~~Figure out what to put in the paper!~~
 * ~~Make overleaf project~~
 * ~~Update David Roberts' slides to fix errors~~
 * ~~write mass as integer+fractional~~
 * ~~Indicate when the list of fields in a family is complete, and handle no fields more gracefully~~
 * ~~Space between brackets in empty list~~
 * ~~decimal mass, mass_missing for sorting and searching~~
 * ~~Write lots of knowls~~
   * ~~Field count~~
 * ~~Update data to divide rams by p-1~~
 * ~~"Abs. Artin slopes" rather than "visible slopes" in search columns~~
 * ~~Talk about what to do with hollow green squares: they are the only things that do not correspond to coefficients in the generic polynomial.~~
 * ~~Should we replace the grid display on the front page with Counts search like in modular forms?~~ No, but we link to such a search.
 * ~~Add more degree 16 extensions of Q2, presumably by finding more Galois splitting models~~
 * ~~Once we have a p-adic polredabs, update the defining polynomials in lf_fields to use it, run in lf_families to get a collection of defining polynomials~~
 * ~~Write first paper~~
 * ~~Finish p-adic polredabs code~~
 * ~~Update grey bands to get darker with more overlap~~
 * ~~Replace "Num. poly" with "Ambiguity," which is the ratio of Num poly by the mass.  It will be a divisor of the degree, and an upper bound for the number of automorphisms for any field in the family.  Equal to p^(num red dots) * f (include base_aut?).  Be able to search on it.~~
 * ~~Change "mass missing" to "mass found"~~
 * ~~Port the 2-d dynamic stats table to the family page itself so you don't have to go to another page; deal with large tables better (like [2.2.2_5_7_9](https://olive.lmfdb.xyz/padicField/dynamic_stats?p=2&n=16&visible_quantifier=exactly&visible=[2,+7%2F2,+9%2F2,+11%2F2]&col1=slopes&totals1=yes&col2=galois_label&totals2=yes&proportions=none).  We could make sub-tables organized by the size of the Galois group.~~
 * ~~For families of large degree, labels on Eisenstein diagram [overlap](https://olive.lmfdb.xyz/padicField/family/2.1.44.130a)~~
 * ~~Add a version of Figure 2.2 from our paper to the family page~~
 * ~~Add Underlying data links~~
 * ~~Total only shows up some subtables for the packet display (e.g. [here](https://olive.lmfdb.xyz/padicField/family/2.1.16.71a); add better border above/left of total row/col; make it clearer how subtables are defined; make row headers be hidden slopes rather than all slopes~~
 * ~~In the "Fields" section of relative families, many of the column headers need to be updated to clarify that they are absolute invariants, not relative.~~
 * ~~Need to add tame information to row headers in the packet section, in the hidden slopes part of the Varying section and in search results.~~
 * ~~Convert `a/b` to `\frac{a}{b}` everywhere, aligned with our discussion in the paper.~~
 * ~~Label overlap in Herbrand plot (2.1.32.140z)~~
 * ~~Change ramification polygon for family page back to tame=slope 0~~
 * ~~Why does 2.1.44.76f have no relative constituents?~~
 * ~~Add tabs and a toggle between the three pictures for a family~~
 * ~~For 2.1.4.6a2.1, the Galois splitting model doesn't link to the number field page since it's not in polredabs form~~
 * ~~Add "Visible Swan slopes" and "Swan slope content" columns to field search results, hidden by default.~~
 * ~~Update [lf.top_slope](https://olive.lmfdb.xyz/knowledge/show/lf.top_slope) to make clear it's referring to Artin slopes. (now called top Artin slope in the knowl, and it refers to a knowl defining Artin slopes)~~
 * ~~Add "Visible Swan slopes" to the "Invariants" section of a [field homepage](https://olive.lmfdb.xyz/padicField/2.1.16.71a1.913) and "Wild Swan slopes" to the "Invariants of the Galois closure" section.~~
 * ~~Knowl for lf.tilts only talks about rams~~
