

## LMFDB front end

 * Add more dynamic columns to dynamic stats, fix links (currently if you click on one of the entries [here](https://olive.lmfdb.xyz/padicField/dynamic_stats?p=2&n=8&visible_quantifier=exactly&visible=[2,+3,+17%2F4]&col1=slopes&totals1=yes&col2=galois_label&totals2=yes&proportions=none) it also includes fields from another family since visible isn't being unparsed correctly).
 * Port the 2-d dynamic stats table to the family page itself so you don't have to go to another page; deal with large tables better (like [2.2.2_5_7_9](https://olive.lmfdb.xyz/padicField/dynamic_stats?p=2&n=16&visible_quantifier=exactly&visible=[2,+7%2F2,+9%2F2,+11%2F2]&col1=slopes&totals1=yes&col2=galois_label&totals2=yes&proportions=none).  We should add tame information to the hidden slopes (and remove the visible slopes from the row headers); ideally we could make sub-tables organized by the size of the Galois group.
 * Write lots of knowls
   * the picture
   * the generic defining polynomial
   * the label of a family (this will probably change from the current Serre-Swan slopes to Artin slopes soon to match the rest of the LMFDB, which will break some of the links in this email...)
   * the rams for a family
   * the heights for a family
   * if we intend to systematically add both Serre-Swan and Artin slopes to other pages, various knowls will need to be updated/created.
   * Num fields
   * Num poly (this should be related to the knowl for the generic defining polynomial, and may want to wait to see if we are able to reduce the set of polynomials considered)
 * Indicate when the list of fields in a family is complete, and handle no fields more gracefully
 * Add ramification polygon to family page
 * Add the canonical subfields to the rows in the list of fields at the bottom of each family page.
 * Update grey bands to get darker with more overlap
 * David mentioned work in number fields for finding Galois fields with root discriminant bounded by a specific value (around 45).  This is connected to the compositum game: trying to find different number fields with very similar ramification so that their compositum has small root discriminant.  For this, it would be helpful to be able to search based on the following partial order on number fields: For every prime we associate two numbers: the top wild slope s (which will be 1 for tame and 0 for unramified) and the lcm t of the tame degrees above that prime.  We say K <= L if s(K) <= s(L) and s(K) | s(L) for all primes p.  Given a number field L we want to find all other K in the database that are less than or equal to L.
 * Make tame part have slope -1 in ramification polygon
 * Name tame fields as Q_q(nth_root(pi)) for an appropriate pi, then update the base method in family.py.
 * Fix __iter__ when base is not Qp.
 * Check that the new hidden column is correctly sorted, update hidden_slopes attribute in family.py to use it.
 * Replace "Num. poly" with "Redundancy," which is the ratio of Num poly by the mass.  It will be a divisor of the degree, and an upper bound for the number of automorphisms for any field in the family.  Equal to p^(num red dots) * base_aut.  Be able to search on it.
 * Add slope multiplicities column, make it possible to search.

## LMFDB data

 * Add more degree 16 extensions of Q2, presumably by finding more Galois splitting models
 * Once we have a p-adic polredabs, update the defining polynomials in lf_fields to use it, run in lf_families to get a collection of defining polynomials
 * Once we have a p-adic polredabs, compute the field labels for each polynomial in families and store in another table.
 * add relative Galois groups to the database (Galois group over the maximal unramified subextension and the maximal tame subextension).

## Other code

 * Finish p-adic polredabs code
 * Given a non-Galois K/Qp, write code to find the slope filtration on Gal(K/Qp) as a sequence of subgroups.  John suggested first finding the slopes (giving the sizes for the groups in the filtration), then trying to use resolvents to determine which options have the correct sized fixed field.
 * Given a non-Galois K/Qp, write code to find the canonical filtration by visible slopes Q_p < K_1 < K_2 < ... < K_j < K, where each extension K_j/K_{j-1} has a single slope associated to it, strictly increasing.  Hopefully this can run in moderate degree (e.g. 1000 over Q2), and may be easier than the general find-subfields problem.

## Theory pursuits

 * Investigate how Galois groups vary in a family (is there a generic group?  can we predict which groups will arise a priori somehow?)

## Paper

 * Figure out what to put in the paper!
 * Make overleaf project

## Other

 * Review our previous reports and pick up threads that we haven't been working on.
 * Update David Roberts' slides to fix errors

## Done

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
