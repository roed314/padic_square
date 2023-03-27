import functools

def ramification_polygon(psi):
        r"""
        Returns the ramification polygon of `psi`.

        The ramification polygon is the Newton polygon of psi(a*x+a)/a^n where a is a root of `psi` and `n` is the degree of `psi`.

        EXAMPLES::

        The vertices of a ramification polygon and the slopes of its segments::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^27+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6
            sage: rp = ramification_polygon(f)
            sage: rp.vertices()
            [(1, 51), (3, 24), (9, 9), (27, 0)]
            sage: rp.slopes(repetition=False)
            [-27/2, -5/2, -1/2]

        A ramification polygon with a horizontal segment::

            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^108+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6
            sage: rp = ramification_polygon(f)
            sage: rp.vertices()
            [(1, 132), (3, 24), (9, 9), (27, 0), (108, 0)]
            sage: rp.slopes(repetition=False)
            [-54, -5/2, -1/2, 0]

        AUTHORS:

        - Brian Sinclair and Sebastian Pauli (2017-07-19): initial version

        """
        #if not psi.is_eisenstein():
        #    raise ValueError("ramification polygons are only defined for Eisenstein polynomials and this polynomial is not Eisenstein")

        from sage.geometry.newton_polygon import NewtonPolygon

        # First we find the ordinates of points above p^k
        verts = []
        vv = [cc.valuation() for cc in psi]
        k = psi.base_ring()
        p = k.prime()
        n = ZZ(psi.degree())
        su = n.valuation(p)
        abscissa = 1
        for i in range(su):
            abscissa = p**i
            ordinate = min([n * (k(binomial(kk,abscissa)).valuation() + vv[kk] - 1) + kk for kk in range(abscissa,n+1)])
            verts.append((abscissa,ordinate))

        # Now we add the tame segment
        for i in range(p**su,n):
            if binomial(n,i).valuation(p) == 0:
                verts.append((i,0))

        # Finally the point for the monic leading term
        verts.append((n,0))

        return NewtonPolygon(verts)

def ramification_polygon_with_colinear_points(psi):
        r"""
        Returns the ramification polygon of `psi` as a list of points including all points on segments of the lower convex hull.

        The ramification polygon is the Newton polygon of `psi(a*x+a)/a^n` where `a` is a root of psi and `n` is the degree of `psi`.

        EXAMPLES::

        We compare the output of this method with the output of `ramification_polygon`.

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(5,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^25 + 20*x^6 + 20*x^5 + 5
            sage: f.ramification_polygon_with_colinear_points()
            [(1, 6), (5, 5), (25, 0)]

        The colinear point (5, 5) is missing on the ramification polygon returned by::

            sage: f.ramification_polygon()
            Finite Newton polygon with 2 vertices: (1, 6), (25, 0)

        When the generated extensions has a tamely ramified subextension::

            sage: f = x^100 + 20*x^6 + 20*x^5 + 5
            sage: f.ramification_polygon_with_colinear_points()
            [(1, 6), (5, 5), (25, 0), (50, 0), (75, 0), (100, 0)]
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 6), (25, 0), (100, 0)

        AUTHORS:

        - Brian Sinclair (2017-07-19): initial version

        """
        #if not psi.is_eisenstein():
        #    raise ValueError("the polynomial psi must be Eisenstein")

        # First we find the ordinates of points above p^k
        verts = []
        vv = [cc.valuation() for cc in psi]
        k = psi.base_ring()
        p = k.prime()
        n = ZZ(psi.degree())
        su = n.valuation(p)
        abscissa = 1
        for i in range(su):
            abscissa = p**i
            ordinate = min([n * (k(binomial(kk,abscissa)).valuation() + vv[kk] - 1) + kk for kk in range(abscissa,n+1)])
            verts.append((abscissa,ordinate))

        # Now we add the tame segment
        for i in range(p**su,n):
            if binomial(n,i).valuation(p) == 0:
                verts.append((i,ZZ(0)))

        # Finally the point for the monic leading term
        verts.append((n,ZZ(0)))

        # Next we need to take the lower convex hull of these points
        def cross(o, a, b):
            """
            2D cross product of the vectors oa and ob.
            """
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        lower = [verts[0]]
        segments = []
        for i in range(1,len(verts)):
            # We check cross < 0 since we want to retain points on the boundary.
            while len(lower) >= 2 and cross(lower[-2], lower[-1], verts[i]) < 0:
                lower.pop()
            lower.append(verts[i])
        if len(lower) <= 1:
            raise ValueError("Not enough vertices")
        return lower

def hasse_herbrand(psi,m):
        r"""
        Returns `n` times the (generalized) Hasse-Herbrand function of `psi` evaluated at `m`.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(2,200,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^16+2
            sage: ramification_polygon(f)
            Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)

        We evaluate the Hasse-Herbrand function at various integers.

            sage: [hasse_herbrand(f,m) for m in range(18)]
            [0, 16, 32, 40, 48, 52, 56, 60, 64, 66, 68, 70, 72, 74, 76, 78, 80, 81]

        Now a different degree:

            sage: f = x^80+2
            sage: ramification_polygon(f)
            Finite Newton polygon with 6 vertices: (1, 320), (2, 240), (4, 160), (8, 80), (16, 0), (80, 0)
            sage: [f.hasse_herbrand(m) for m in range(16)]
            [0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 168, 176, 184, 192, 200]

        AUTHORS:

        - Brian Sinclair (2017-07-20): initial version
        """
        return min([pt[1]+m*pt[0] for pt in ramification_polygon(psi).vertices()])

def residual_polynomials(psi):
        r"""
        Returns a list of the residual polynomials of the ramification polynomial of an Eisenstein polynomial psi.

        EXAMPLES::

        The residual polynomials of the segments of the ramification polygon of an Eisenstein polynomial::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9+6*x^3+3
            sage: f.is_eisenstein()
            True
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
            sage: f.residual_polynomials()
            [z + 2, z^3 + 1]

        A ramfication polygon with a horizontal segment::

            sage: f = x^90+6*x^3+3
            sage: f.ramification_polygon()
            Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-45, -1/2, 0]
            sage: f.residual_polynomials()
            [z^2 + 2, z^3 + 1, z^81 + z^72 + 1]

        A ramification polygon with more segments::

            sage: R = ZpFM(2,300,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^16+2
            sage: f.ramification_polynomial()
            (1 + O(a^4800))*y^16 + (16 + O(a^4800))*y^15 + (120 + O(a^4800))*y^14 + (560 + O(a^4800))*y^13 + (1820 + O(a^4800))*y^12 + (4368 + O(a^4800))*y^11 + (8008 + O(a^4800))*y^10 + (11440 + O(a^4800))*y^9 + (12870 + O(a^4800))*y^8 + (11440 + O(a^4800))*y^7 + (8008 + O(a^4800))*y^6 + (4368 + O(a^4800))*y^5 + (1820 + O(a^4800))*y^4 + (560 + O(a^4800))*y^3 + (120 + O(a^4800))*y^2 + (16 + O(a^4800))*y + (0 + O(a^4800))
            sage: f.ramification_polygon()
            Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)
            sage: f.residual_polynomials()
            [z + 1, z^2 + 1, z^4 + 1, z^8 + 1]

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        """
        #if not psi.is_eisenstein():
        #    raise ValueError("the polynomial psi must be Eisenstein")

        n = psi.degree()
        phi0 = psi.constant_coefficient()
        phi00 = phi0>>1
        Rx = psi.parent()
        R = Rx.base()
        p = R.prime()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()

        rp = ramification_polygon_with_colinear_points(psi)
        respols = []
        j = 0
        while j<len(rp)-1:
            k = j
            p_s_k = rp[k][0]
            slope = (rp[j+1][1]-rp[j][1])/(rp[j+1][0]-rp[j][0])
            e = slope.denominator()
            thispol = Fz(0)
            while True:
                p_s_i = rp[j][0]
                a_i, b_i = rp[j][1].quo_rem(n)
                if b_i == 0:
                    a_i -= 1
                    b_i  = n
                thispol += ((psi[b_i]*binomial(b_i,p_s_i)*(-phi00)**(-a_i-1))>>(a_i+1)).residue()*z**((p_s_i-p_s_k)//e)
                if j>=len(rp)-1 or (rp[j+1][1]-rp[j][1])/(rp[j+1][0]-rp[j][0]) != slope:
                    break
                j+=1
            respols.append(thispol)

        return respols

def residual_polynomial_of_component(psi,m):
        r"""
        Return the residual polynomial ``S_m`` of the (`-m`)-component of the ramifation polygon of polynomials psi, which must be Eisenstein.

        Let N be the ramification polygon then {(k,w) in N | (`-m`)k + w = min{(`-m`)l+u|(l,u) in N} } is the (`-m`)-component of N.

        INPUT::

            A natural number `m`

        OUTPUT::

            The residual polynomial of the (`-m`)-component of the ramificaton polygon of psi.

        EXAMPLES::

        In our first example, we have a polynomial over the degree 2 unramified extension of ``\QQ_2``
        which has a ramification polygon with two segments of integral slope::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R.<g> = ZqFM(4,30); Rx.<x> = R[]
            sage: f = x^8 + 2*g*x^6 + 4*g*x + 2
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 9), (2, 6), (8, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-3, -1]
            sage: f.residual_polynomials()
            [g0*z + g0, z^6 + g0]
            sage: [f.residual_polynomial_of_component(m) for m in range(1,10)]
            [z^8 + g0*z^2, z^2, g0*z^2 + g0*z, z, z, z, z, z, z]

        Here we have a nonic polynomial over ``\QQ_3`` whose ramification polygon has no segemnts
        of integral slope::

            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9 + 6*x^3 + 3
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-9/2, -1/2]
            sage: f.residual_polynomials()
            [z + 2, z^3 + 1]
            sage: [f.residual_polynomial_of_component(m) for m in range(0,10)]
            [z^9, z^3, z^3, z^3, z^3, z, z, z, z, z]

        Here the ramification polygon has a horizontal segment::

            sage: f = x^90 + 6*x^3 + 3
            sage: f.ramification_polygon()
            Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-45, -1/2, 0]
            sage: f.residual_polynomials()
            [z^2 + 2, z^3 + 1, z^81 + z^72 + 1]
            sage: [f.residual_polynomial_of_component(m) for m in range(0,10)]
            [z^90 + z^81 + z^9, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3]

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        REFERENCES:

        [PS17] S. Pauli and B. Sinclair, "Enumerating Extensions of (pi)-adic Fields with Given Invariants", International Journal of Number Theory
          (2017)
        """

        #if not psi.is_eisenstein():
        #    raise ValueError("residual polynomials are only defined for Eisenstein polynomials")

        Rx = psi.parent()
        R = Rx.base()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()

        rp = ramification_polygon(psi)

        if -m in rp.slopes(repetition=False):
            i = rp.slopes(repetition=False).index(-m)
            return residual_polynomials(psi)[i]*z**rp.vertices()[i][0]
        else:
            L = [v[1]+v[0]*m for v in rp.vertices()]
            mini = min(L)
            mindex = L.index(mini)
            return z**rp.vertices()[mindex][0]

def deformed_eisenstein(psi, m, theta, trunc, test=False):
        r"""
        Return the "deformed" Eisenstein polynomial, ie. the minimal polynomial of the uniformizer `q` such that `q + theta*q^(m+1) = pi`.
        Note that we also have that `q = pi - theta*pi^(m+1) + O(pi^(m+2))`

        EXAMPLES::

        We deform an Eisenstein polynomial and check the relationship between `q` and `pi`.

            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^9+9*x^2+3
            sage: m = 2; theta = 1; trunc = 180
            sage: g = deformed_eisenstein(f,m,theta,trunc)
            sage: g
            sage: S.<q> = R.ext(g)
            sage: pi = q + theta*q^(m+1)
            sage: h = f.change_ring(S)
            sage: h(pi).is_zero()

        AUTHORS:

        - Maurizio Monge (2014-11-14): initial version
        - Sebastian Pauli and Brian Sinclair (2017-07-20): bang to sage 8.0, documentation

        """
        if theta==0:
           return psi

        #print("deform",psi,"by",theta,"* pi ^",m)
        #if not psi.is_eisenstein():
        #    raise ValueError("this function can only deform Eisenstein polynomials")
        #print("deforming")
        f = psi
        n = f.degree()
        x = f.variables()[0]
    
        OKx = f.parent()
        OKu.<u> = PolynomialRing(OKx)
        fu = f(u)
        gu = u-theta*u**(m+1)
        ru = fu.resultant(x-gu)
        #print("resultant",ru)
        return ru

        # truncation level, because x^trunc = 0 mod piK^prec, as v(x) = 1/n
        g = f(x + theta*x**(m+1)).truncate(trunc) # IMPROVE-ME: use f_0 instead of x^n
        if test:        
          if not generate_isomorphic(psi,g): 
             raise ValueError("transformation failed")
        #print("g",g)
        prev = (0, n)
        while(g.degree() > n):
            # where f has terms we want to kick out
            extra_range = range(n+1, g.degree()+1)

            cur_prec = min([(g[i].valuation())*n + i for i in extra_range])
            if cur_prec >= trunc:
                break

            # lexicographically minimal pair (val(f_i), i), for i in the extra range
            min_val, min_idx = min([(g[i].valuation(), i) for i in extra_range])

            # the bad terms should be going away...
            assert(prev < (min_val, min_idx))
            prev = (min_val, min_idx)

            # subtract, from g, (badmononial/x^n)*g
            g = (g - g[min_idx] * g.shift(min_idx-n)).truncate(trunc)
        
            if test:        
             if not generate_isomorphic(psi,g): 
                raise ValueError("deformation loop failed at n =",n)

        gdef = g.truncate(n+1).monic()
        
        if test:        
          if not generate_isomorphic(psi,g): 
             raise ValueError("deformation failed")

        return gdef 





def monge_reduce(psi,all=False,test=False):
        r"""
        Return the Monge-reduced polynomial that generates an extensions isomorphic to the extensions generated by the Eisenstein polynomial `psi`.

        When fixing set of representatives for the classes of elements of the residue class field occuring in the Monge reduction,
        the Monge-reduced polynomials are unique.  We make the following choices.

        If the coefficient ring of `psi` is unramified, we choose the representatives of the classes of elements of the residue class field 
        such that their absolute value is minimal.  Otherwise we order the representatives lexicographically.
        Note that this representation depends on the generating polynomial of the unramified part of the extension,

        If all == True then all possible Monge reductions of the polynomial are returned.

        EXAMPLES::

        We Monge-reduce a polynomial,.

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9+249*x^3+486*x+30
            sage: g = monge_reduce(f,test=True)
            sage: g
            (1 + O(3^30))*x^9 + ... + (2*3 + O(3^30))*x^3 + ... + (3 + O(3^30))

        If the polynomial is Monge-reduced it does not change when reduced again::

            sage: h = monge_reduce(f,test=True)
            sage: h == g
            True

        The Monge-reduction of a polynomial generating a tamely ramified extension::

            sage: f = x^20+249*x^3+486*x+30
            sage: g = monge_reduce(f,test=True)
            sage: g
            (1 + O(3^30))*x^20 + ... + (3 + O(3^30))

        The Monge-reduction of a polynomial generating a tamely ramified extension of large degree::

            sage: f = x^90+249*x^81+486*x^18+30
            sage: g = monge_reduce(f,test=True)
            sage: ZZ['X'](g)
            X^90 + 3*X^81 + 9*X^78 + 9*X^72 + 9*X^54 + 54*X^44 + 54*X^43 + 54*X^41 + 54*X^40 + 18*X^39 + 54*X^38 + 27*X^26 + 54*X^13 + 3

        We use Monge reduction to verify that two polynomials generate isomorphic extensions

            sage: R = ZpFM(5,20); Rx.<x> = R[]
            sage: f = x^25+15625*x^4+5
            sage: g = x^25+5
            sage: monge_reduce(f,test=True) == monge_reduce(g)
            True

        Monge-reduction over an unramified extensions::

            sage: R.<g> = ZqFM(4,30); Rx.<x> = R[]
            sage: f = x^8 + 66*g*x^6 + 132*g*x + 258
            sage: monge_reduce(f,test=True)
            (1 + O(2^30))*x^8 + (O(2^30))*x^7 + (g*2 + O(2^30))*x^6 + ... + (g*2^2 + O(2^30))*x + (2 + O(2^30))

        Monge reduction over a totally ramified extension::

            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: S.<a> = R.ext(x^3+9*x+3); Sy.<y> = S[]
            sage: f = y^6+6*y^2+a
            sage: monge_reduce(f)
            (1 + O(a^90))*y^6 + (O(a^90))*y^5 + (O(a^90))*y^4 + (O(a^90))*y^3 + (a^3 + O(a^90))*y^2 + (O(a^90))*y + (a + O(a^90))

            sage: R.<t> = ZqFM(4,500); Rx.<x> = R[]
            sage: f = x^8 + 4*x^7 + 12*x^6 + 20*x^5 + (-2*t + 22)*x^4 + (-8*t + 16)*x^3 + (-8*t + 8)*x^2 - 8*t*x - 2*t
            sage: monge_reduce(f, test=True)

            sage: R = ZpFM(2,30);Rx.<x> = R[]
            sage: f = x^2+14
            sage: monge_reduce(f,test=True)


        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        REFERENCES:

        [Mon14] M. Monge, "A family of Eisenstein polynomials
          generating totally ramified extensions, identification of extensions and
          construction of class fields." International Journal of Number Theory
          (2014): 1-29.
        """

        #if not psi.is_eisenstein():
        #    raise ValueError("only Eisenstein polynomials can be Monge-reduced")

        f = psi
        n = f.degree()
        J = f.discriminant().valuation()
        RT = f.parent()
        T = RT.gen()
        R = RT.base()
        p = R.prime()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()
        r = F.degree()
        xi = F.primitive_element() 
        q = F.order()

        def canonical_representative_mult(alpha,k):
            ai = F(alpha)**(q-2)
            if F.is_prime_field():
                bset = [ZZ(alpha*s) for s in [a**k for a in F if not a.is_zero()]]
                return F(min(bset))
            else:
                for i in range(0,q-1):
                    #print("xi^i",xi^i)
                    h = z**k-ai*xi**i
                    if len(h.roots())>0:
                       return xi**i

        def canonical_representative_add(alpha, funct):
            #print("canonical_representative_add")
            if F.is_prime_field():
                #print("prime field")
                bset = [ZZ(alpha+s) for s in [funct(a) for a in F]]
                return F(min(bset)) 
            for i in range(0,q-1):
                #print("extension")
                h = alpha+funct-xi**i
                #print("xi^i",xi^i)
                if len(h.roots())>0:
                    return xi**i

        def solve_naive_obsolete(funct,gamma):
        # solve poly(t)=gamma
        # terrible
            sol = []
            for a in F:
                if funct(a)==gamma:
                   sol.append(a)
                   # return a
            if sol == []:
              raise ValueError("No solution found")
            return sol



        def solve_naive(funct,gamma):
            #print("solving",funct,"=",gamma)
            h = funct-gamma
            r = h.roots()
            #print("solutions",[a[0] for a in r])
            return [a[0] for a in r]


        phi0 = f.constant_coefficient()
        eta = (phi0 >> 1).residue()
        alpha = eta

        #print(" reduction step 0 -- taking care of the constant coefficient"
        #print("n",n,type(n),"p",p,type(p))
        n = ZZ(n)
        if n == p**n.valuation(p):
            beta = 1
        else:
            beta = canonical_representative_mult(eta,n) #set([a**n for a in F if not a.is_zero()]))
        #print("beta",beta)

        #print(" reduction step 0 -- solving")
        #print(set(solve_naive(z**n,beta/alpha))==set(solve_naive_obsolete(z**n,beta/alpha)))
        Thetas = solve_naive(z**n, beta/alpha)
        #print("theta",Thetas)
        # Theta = R(theta)
        # f = theta**n*f(theta^-1*T)
         
        # list of Monge reduced polynomials
        # M = {RT([f[i]*R(theta)**(n-i) for i in range(0,n+1)]) for theta in Thetas}
        M = [RT([f[i]*R(theta)**(n-i) for i in range(0,n+1)]) for theta in Thetas]
        #print("M",M)

        if test:
          for psi in M:
            if not generate_isomorphic(psi,f): 
               raise ValueError("constant coefficient reduction failed")
        #print("other reduction steps")

        def f_ij(f,m):
            lev = hasse_herbrand(psi,m)
            i = lev % n
            j = ZZ((n - i + lev) // n)
            fij = F(f[i].expansion(j))
            return fij, i, j


        J0 = ramification_polygon(psi).vertices()[0][1]
        trunc = n + 2*J0 # truncate here in deform
        for m in range(1,(-min(ramification_polygon(psi).slopes())).ceil()):
            #print("m",m)
            alpha, i, j = f_ij(f,m)
            Sm = residual_polynomial_of_component(f,m)
            #beta = canonical_representative_add(alpha,[eta**j*Sm(a) for a in F])
            beta = canonical_representative_add(alpha,eta**j*Sm)
            #print("fbeta",beta)
            #print("Thetas: solutions of",eta**j*Sm," = ",alpha,"-",beta)
            Thetas = solve_naive(eta**j*Sm,alpha-beta)
            #print(set(solve_naive(eta**j*Sm,alpha-beta))==set(solve_naive_obsolete(eta**j*Sm,alpha-beta)))
            # M = {deformed_eisenstein(ff,m, -R(theta), trunc,test=test) for theta in Thetas for ff in M}
            M = [deformed_eisenstein(ff,m, -R(theta), trunc,test=test) for theta in Thetas for ff in M]
            #print(m,"M",M)

            # f = deformed_eisenstein(f,m, -R(theta), trunc,test=test)

        #print("M",M)
        #print("f",f)
        if test:
          #print("before beyond")        
          for psi in M:
            if not generate_isomorphic(psi,f): 
               raise ValueError("Sm reduction failed")
        
        # Find the last break in the Hasse-Herbrand function of psi
        hhslope = n
        m = 1
        hhm = 0
        while hhslope > 1:
            #print("slope",hhslope)
            lrb = hhm
            hhm = hasse_herbrand(psi,m)
            hhslope = hhm - lrb
            m += 1

        # If f generates a tamely ramified extension, then lrb == 0, which
        # would clear all coefficients.  Set to 1, it preserves f_{0,1}.
        lrb = max(lrb,1)
        
        # All coefficients lexicographically beyond
        #    (lrb mod n, (n-(lrb mod n) + lrb)/n)
        # can be set to zero.
        
        i = lrb % n
        j = ZZ((n - i + lrb) // n)
        #print("Krasner bound reduction")
        #print("kM",M)
        
        if R.ramification_index() == 1:
            #print("e=1")
            # f = RT([cc.add_bigoh(j+1) for cc in f.list()[:i+1]] + [cc.add_bigoh(j) for cc in f.list()[i+1:n]] + [R.one()])
            # M = {RT([cc.add_bigoh(j+1) for cc in ff.list()[:i+1]] + [cc.add_bigoh(j) for cc in ff.list()[i+1:n]] + [R.one()]) for ff in M}
            M = [RT([cc.add_bigoh(j+1) for cc in ff.list()[:i+1]] + [cc.add_bigoh(j) for cc in ff.list()[i+1:n]] + [R.one()]) for ff in M]
        else:
            # M = {RT([cc - (cc >> j+1 << j+1) for cc in ff.list()[:i+1]] + [cc - (cc >> j << j) for cc in ff.list()[i+1:n]] + [R.one()]) for ff in M}
            M = [RT([cc - (cc >> j+1 << j+1) for cc in ff.list()[:i+1]] + [cc - (cc >> j << j) for cc in ff.list()[i+1:n]] + [R.one()]) for ff in M]
                
        #print("M",M)
        if test:        
          #print("after beyond")        
          for psi in M:
            if not generate_isomorphic(psi,f): 
               raise ValueError("krasner reduction failed")
        
        # return all monge reduction of f

        if all: 
          return M

        # return the nicest one:

        ZX = PolynomialRing(Integers(),'X')
        nz = min([ff.number_of_terms() for ff in M])
        # M = [ZX(ff) for ff in M if ff.number_of_terms() == nz] 
        M = [ff for ff in M if ff.number_of_terms() == nz] 
        M.sort()
        # M = [RT(ff) for ff in M]
        return M[0]

def content(psi):
        r"""
        Compute the content of this polynomial.

        OUTPUT:

        If this is the zero polynomial, return the constant coefficient.
        Otherwise, since the content is only defined up to a unit, return the
        content as `\pi^k` with maximal precision where `k` is the minimal
        valuation of any of the coefficients.

        EXAMPLES::

            sage: K = Zp(13,7)
            sage: R.<t> = K[]
            sage: f = 13^7*t^3 + K(169,4)*t - 13^4
            sage: f.content()
            13^2 + O(13^9)
            sage: R(0).content()
            0
            sage: f = R(K(0,3)); f
            (O(13^3))
            sage: f.content()
            O(13^3)

            sage: P.<x> = ZZ[]
            sage: f = x + 2
            sage: f.content()
            1
            sage: fp = f.change_ring(pAdicRing(2, 10))
            sage: fp
            (1 + O(2^10))*x + (2 + O(2^11))
            sage: fp.content()
            1 + O(2^10)
            sage: (2*fp).content()
            2 + O(2^11)

        Over a field it would be sufficient to return only zero or one, as the
        content is only defined up to multiplication with a unit. However, we
        return `\pi^k` where `k` is the minimal valuation of any coefficient::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: f = 13^7*t^3 + K(169,4)*t - 13^-4
            sage: f.content()
            13^-4 + O(13^3)
            sage: f = R.zero()
            sage: f.content()
            0
            sage: f = R(K(0,3))
            sage: f.content()
            O(13^3)
            sage: f = 13*t^3 + K(0,1)*t
            sage: f.content()
            13 + O(13^8)

        """
        if psi.is_zero():
            return psi[0]
        else:
            try:
               val = psi.valuation(val_of_var=0)
            except TypeError:
               val = min([c.valuation() for c in psi])
            return psi.base_ring()(psi.base_ring().uniformizer_pow(val))


def normalize(psi,monic_if_possible=False):
        r"""
        Return a polynomial and divisor such that polynomial = psi // divisor and
        the polynomial is such that 0 is the lowest valuation amongst the coefficients.

        INPUT::

        - ``monic_if_possible`` Use the leading coefficient as the scaling divisor if possible.

        OUTPUT::

            Normalized polynomial and the scaling divisor

        EXAMPLES::
            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30,print_mode='terse'); Rx.<x> = R[]
            sage: f = 9*x^9+3*x^2+27
            sage: f.normalized()
            ((3 + O(3^30))*x^9 + ... + (9 + O(3^30)), 3 + O(3^30))

        Here we see the unit factored out of the leading term if possible::

            sage: g = 6*x^9+3*x^2+27
            sage: g.normalized()
            ((2 + O(3^30))*x^9 + ... + (1 + O(3^30))*x^2 + (0 + O(3^30))*x + (9 + O(3^30)), 3 + O(3^30))
            sage: g.normalized(monic_if_possible=True)
            ((1 + O(3^30))*x^9 + ... + (102945566047325 + O(3^30))*x^2 + (0 + O(3^30))*x + (102945566047329 + O(3^30)), 6 + O(3^30))

        """
        normalize_by = content(psi)
        psi_normal = psi.parent()([c>>normalize_by.valuation() for c in psi])
        if monic_if_possible and psi_normal.leading_coefficient().valuation() == 0:
            # if valuation of leading coefficient of normalized polynomial will be 0
            # normalize by leading coefficient
            normalize_by = normalize_by*psi_normal.leading_coefficient()
            psi_normal = psi.parent()(psi_normal/psi_normal.leading_coefficient())
        return psi_normal,normalize_by

def has_root(phi):
 
       f = phi
       n = f.degree()
       Opx = f.parent()
       x = Opx.gen()
       Op = Opx.base()
       pi = Op.uniformizer()
       Fq = Op.residue_class_field()
       Fqy = PolynomialRing(Fq,"y")
       y = Fqy.gen()

       psi = normalize(phi)[0] 
       #print("psi",psi)
       C = [ (psi,Op(0),0) ]  

       while len(C) !=0:
         #print("C",C)
         this = C.pop()
         #print("this",this)
         (psi,delta,s) = this
         #print("s",s)
         #print("delta",delta)
         #print("psi(delta)",psi(delta))
         psib = Fqy(psi)
         betas = psib.roots()
         #print([a[1] for a in betas])
         #print("betas",betas)
         for beta in betas:
           #print("delta",delta)
           #print("beta",beta)
           newpsi = psi(pi*x+Op(beta[0]))
           newpsi = normalize(newpsi)[0]
           #print(s,newpsi.degree())
           if beta[1] == 1:
              #print("single root",betas)
              return True, delta+Op(beta[0])*pi**s
           #if s > pi.precision_absolute():
           #  R.append(delta + Op(beta[0])*pi**s)
           elif newpsi.degree() == 1:
              #print("degree 1")
              return True, delta+Op(beta[0])*pi**s
           #  gamma = newton_lift(newpsi,Op(beta[0]))
           #  R.append(delta + gamma*pi**s)
           else:
             C.append((newpsi,delta+Op(beta[0])*pi**s,s+1))
       return False, None
      
def generate_isomorphic(phi,psi):
           #print("testing")
           RT = psi.parent()
           R = RT.base()
           S.<a> = R.ext(phi)
           St.<t> = PolynomialRing(S)
           return has_root(St(psi))[0]

