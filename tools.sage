# composition m2 * m1
def compose(m1, m2):
  return (m2[0].subs(x=m1[0],y=m1[1]), m2[1].subs(x=m1[0],y=m1[1]))

# endomorphism [r] on E. nr = norm_{K|Q}(r), hpoly: the Hilbert class polynomial for E
def mult_r(E, r, nr, hpoly):
    factors = factor(nr)
    idxs = None
    for f in factors:
        if f[1] > 1:
            tmp = cartesian_product([[0..1] for _ in range(f[1])])
        else:
            tmp = [0..1]
        if idxs == None:
            idxs = tmp
        else:
            idxs = cartesian_product([idxs, tmp])

    for idx in idxs:
        isogs = []
        Ed = E
        i = 0
        for l, e in factors:
            for _ in range(e):
                horizontals = [phi for phi in Ed.isogenies_prime_degree(l) if hpoly(phi.codomain().j_invariant()) == 0]
                isog = horizontals[idx[i]]
                Ed = isog.codomain()
                isogs.append(isog)
                i += 1

        if Ed.j_invariant() == E.j_invariant():
            isom = Ed.isomorphism_to(E)
            isogs[-1].set_post_isomorphism(isom)
            alpha = isogs[0].rational_maps()
            for isog in isogs[1:]:
                alpha = compose(alpha, isog.rational_maps())
            diff = coeff_invariant_differential(alpha)
            if diff == r:
                return alpha
            elif diff == - r:
                return (alpha[0], -alpha[1])

# return alpha*w / w, where w = dx/(2y) is the invariant differential
def coeff_invariant_differential(alpha):
    a_x = alpha[0](x=T, y=1)
    a_y = alpha[1](x=T, y=1)
    return a_x.derivative()/a_y

# quadratic residue symbol (a|pk) by the product formula of Hilbert symbols
# (a|pk) = prod_{t|(2)} (a,pk|t) * prod_{I_ell|(a), not I_ell|2}(pk|I_ell),
# where t and I_ell's are prime ideals.
def quadratic_residue_H(a, pk):
    s = 1
    for t in H.primes_above(2):
        s *= H.hilbert_symbol(a, pk, t)
    for I_ell, ex in factor(H.ideal(a)):
        if I_ell.reduce(pk) == 0:
            s = 0
            break
        if I_ell.divides(2):
            continue
        s *= I_ell.residue_symbol(pk, 2)^ex
    return s

# the period of a sequence (a|pk) for k.
def qr_period_k(a):
    period = 2 # The hilbert symbol (a,pk|t) for a prime ideal t|(2) is determined by k mod 2.
    for I_ell, ex in factor(H.ideal(a)):
        if (not I_ell.divides(2)):
            period = lcm(period, norm(I_ell) - 1)
    return period

# return whether list L has period k
def check_period(L, k):
    n = len(L)
    for i in range(n-k):
        if not L[i] == L[i+k]:
            return False
    return True

# return the shortest period of L
def reduce_period(L, period):
    ret = period
    for f in factor(period):
        l, ex = f
        for i in range(ex):
            if check_period(L, ret//l^(ex-i)):
                ret //= l^(ex-i)
                break
    return ret

# square root of t mod p. We assume p is a prime not 1 mod 8.
def sqrt_mod(t, p):
    assert not (p % 8) == 1

    if (p % 4) == 3:
        return power_mod(t, (p + 1)//4, p)
    else:
        s = power_mod(t, (p + 3)//8, p)
        if (s^2 - t) % p == 0:
            return s
        else:
            return s*power_mod(2, (Fk - 1)//4, p)

