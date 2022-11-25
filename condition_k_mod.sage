# the condition that pk is prime to disc(E)
def calcT0(E):
    disc = E.discriminant()
    Idisc = H.ideal(disc)
    period = qr_period_k(disc)

    values = []
    for k in range(1, period+1):
        pk = p_k_mod(k, Idisc).lift()
        values.append(Idisc.is_coprime(pk))
    period = reduce_period(values, period)

    ret = []
    for k in range(period):
        if values[k] == 1:
            ret.append(int((k+1) % period))

    return {"period":int(period), "values":ret}

# norm of H over K
def normHK(lam):
    return to_rel(lam).relative_norm()

# epsilon_tau
def epsilon(lam):
    rd = to_rel(r)
    if K.ideal(4).reduce(lam^3 - 1) == 0 or K.ideal(4).reduce(lam^3 - (-2*rd + 1)) == 0:
        return 1
    elif K.ideal(4).reduce(lam^3 + 1) == 0 or K.ideal(4).reduce(lam^3 + (-2*rd + 1)) == 0:
        return -1
    else:
        # error!
        return 0

# return 1 if the Frobenius on Ed is [c1*r + c2], -1 if -[c1*r + c2]
def sign_frobenius(Ed, c1, c2, rX, rY):
    posi = True
    nega = True
    for P in Ed.points():
        if not P[2] == 0:
            if rX.denominator()(X=P[0], Y=P[1]) == 0:
                Q = Ed([0,1,0])
            else:
                Q = Ed([rX(X=P[0], Y=P[1]), rY(X=P[0], Y=P[1])])
            posi = posi and (P == c1*Q + c2*P)
            nega = nega and (P == -c1*Q - c2*P)
    if posi:
        return 1
    elif nega:
        return -1
    else:
        # error!
        return 0

# return (6*gamma|pH)*epsilon_tau(pK)[pK] == Frobenius
def check_gamma(E, gamma, r, c1, c2, end_r_E):
    pK = c1*r + c2
    pH = H.prime_above(pK)

    # the reduction of E mod pH
    Op = H.residue_field(pH)
    Ed = E.reduction(pH)

    # endomorphism [r] on Ed
    OpXY.<X, Y> = PolynomialRing(Op, 2)
    OpXY = FractionField(OpXY)
    rX = OpXY(end_r_E[0])
    rY = OpXY(end_r_E[1])

    # sign of the Frobenius
    s1 = sign_frobenius(Ed, c1, c2, rX, rY)
    # sign of the Frobenius
    e = epsilon(to_rel(pK))
    # sign of the Frobenius
    s2 = pH.residue_symbol(6*gamma, 2) * e

    return s1 == s2

# Weber function gamma3 of E
def calc_gamma3(E, cs, D):
    end_r_E = mult_r(E, r, K(to_rel(r)).norm(), hilbert_class_polynomial(D))

    gamma = (T^2 - (E.j_invariant() - 1728)).roots(multiplicities=False)[0]
    c1, c2 = cs[0]
    if not check_gamma(E, gamma, r, c1, c2, end_r_E):
        gamma = -gamma

    # just verification
    for c1, c2 in cs[1:]:
        assert check_gamma(E, gamma, r, c1, c2, end_r_E)

    return gamma

# the condition that [norm_{H|K}(pk)] is the Frobenius on E mod pk
def calcT1(E, cs, D):
    gamma = calc_gamma3(E, cs, D)

    # exclude the factor 6 from period because (3|pk) is determined by pk mod 3 so k mod 2
    period = qr_period_k(gamma)

    values = []
    for k in range(1, period+1):
        pk = p_k_mod(k, H.ideal(4*6*gamma)).lift() # Hilbert symbol (a,b|t) for t|2 is determined modulo 2^3
        v = epsilon(normHK(pk))
        v *= quadratic_residue_H(6*gamma, pk)
        values.append(v)
    period = reduce_period(values, period)

    ret = []
    for k in range(period):
        if values[k] == 1:
            ret.append(int((k+1) % period))

    return {"period":int(period), "values":ret}

# the generator of E[2, alpha]
def gen_kernel2(E, alpha):
    for P in E.torsion_points():
        if P.order() == 2:
            if alpha[0].denominator()(x=P[0]/P[2], y=0) == 0:
                return P

# the discriminant of the numerator of x-coordinate of the dual isogeny of the isogeny with kernel <P>
def disc_dual_isog(E, P):
    phi_hat = E.isogeny(P)
    phi = phi_hat.dual()
    num = phi.rational_maps()[0].numerator()
    return num(x=T).discriminant()

# the condition that (0, *) in E is not in (2, r)(E mod pk) 
def calcT2(E, r, D):
    alpha = mult_r(E, r, norm(K(to_rel(r))), hilbert_class_polynomial(D))
    P = gen_kernel2(E, alpha)
    disc = disc_dual_isog(E, P)
    period = qr_period_k(disc)

    values = []
    for k in range(1, period+1):
        pk = p_k_mod(k, H.ideal(8*disc)).lift() # Hilbert symbol (a,b|t) for t|2 is determined modulo 2^3
        values.append(quadratic_residue_H(disc, pk))

    period = reduce_period(values, period)

    ret = []
    for k in range(period):
        if values[k] == -1:
            ret.append(int((k+1) % period))

    return {"period":int(period), "values":ret}

# the intersection of T0, T1 and T2
def integrate_results(Ts):
    ret = []
    period = lcm([t["period"] for t in Ts])
    for k in range(period):
        if prod([k % t["period"] in t["values"] for t in Ts]):
            ret.append(k)
    return {"period":int(period), "values":ret}
