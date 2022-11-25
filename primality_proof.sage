# the x-coordinate of 2*(X:Y:Z) on y^2 = x^3 + Ax + B
def ell_dbl(X, Z, A, B):
    retX = X^4 - 2*A*X^2*Z^2 - 8*B*X*Z^3 + A^2*Z^4
    retZ = 4*Z*(X^3 + A*X*Z^2 + B*Z^3)
    return [retX, retZ]

# the x-coordinate of (X1:Y1:Z1) + (X2:Y2:Z2) on y^2 = x^3 + Ax + B. (X3:Y3:Z3) = (X1:Y1:Z1) - (X2:Y2:Z2)
def ell_diffadd(X1, Z1, X2, Z2, X3, Z3, A, B):
    if not X3 == 0:
        retX = Z3*((X1*X2)^2 - 2*A*(X1*X2*Z1*Z2) - 4*B*Z1*Z2*(X1*Z2 + X2*Z1) + (A*Z1*Z2)^2)
        retZ = X3*(X1*Z2 - X2*Z1)^2
    else:
        retX = 2*(X1*X2*(X1*Z2 + X2*Z1) + A*Z1*Z2*(X1*Z2 + X2*Z1) + 2*B*(Z1*Z2)^2)
        retZ = (X1*Z2 - X2*Z1)^2
    return [retX, retZ]

# the x-coordinate of n*(X:Y:Z) on y^2 = x^3 + Ax + B
def ell_ladder(n, X, Z, A, B):
    X0, Z0 = 1, 0
    X1, Z1 = X, Z
    for k in bin(n)[2:]:
        if k == '0':
            X1, Z1 = ell_diffadd(X0, Z0, X1, Z1, X, Z, A, B)
            X0, Z0 = ell_dbl(X0, Z0, A, B)
        elif k == '1':
            X0, Z0 = ell_diffadd(X0, Z0, X1, Z1, X, Z, A, B)
            X1, Z1 = ell_dbl(X1, Z1, A, B)
    return X0, Z0


"""
return Fk is prime or not.
E: elliptic curve with CM by Q(D),
a_mod_Fk(r_Fk, k) returns a mod Fk, where a is a generator of the Hilbert class field of Q(D),
minpoly_a is the minimal polynomial of a over Q.
"""
def primality_proof(E, D, k, Fk, a_mod_Fk, minpoly_a):
    Z_Fk = quotient(ZZ, Fk)

    sqrtD_Fk = Z_Fk(sqrt_mod(D, Fk))
    if not sqrtD_Fk^2 - D == 0:
        return False

    r_Fk = (1 + sqrtD_Fk) * ((Fk + 1)//2) # (1 + sqrt(D))/2 mod Fk
    a_Fk = a_mod_Fk(r_Fk, k) # a mod Fk

    if not minpoly_a(a_Fk) == 0:
        sqrtD_Fk = -sqrtD_Fk
        r_Fk = (1 + sqrtD_Fk) * ((Fk + 1)//2) # (1 + sqrt(D))/2 mod Fk
        a_Fk = a_mod_Fk(r_Fk, k) # a mod Fk
    assert minpoly_a(a_Fk) == 0

    ar_Fk = a_Fk - r_Fk

    A = Z_Fk(E.a4().polynomial()(x = ZZ(ar_Fk)))
    B = Z_Fk(E.a6().polynomial()(x = ZZ(ar_Fk)))

    P = ell_ladder(c_k(k), 0, 1, A, B)
    for _ in range(6*k - 1):
        P = ell_dbl(P[0], P[1], A, B)
        if gcd(P[1], Fk) > 1:
            return False

    return ell_dbl(P[0], P[1], A, B)[1] == 0
