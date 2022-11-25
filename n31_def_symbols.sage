"""
    Number fields
    K: Q(sqrt{-31})
    H: the Hilbert class field of Q(sqrt{-31})
    r = (1 + sqrt{-31})/2
    a: a root of T^3 + T + 1
    ar = a - r, a primitive element of H/Q
    Hrel: tower of extension H/K/Q
    ar_rel, r_rel = ar, r in Hrel
"""
H.<ar> = NumberField(x^6 + 3*x^5 + 29*x^4 + 55*x^3 + 223*x^2 + 151*x + 379)
HT.<T> = PolynomialRing(H)
for rt in (T^2 - T + 8).roots(multiplicities=False):
    if (ar + rt).minpoly().degree() == 3:
        a = ar + rt
        r = rt
assert(a^3 + a + 1 == 0)
Hrel.<ar_rel,r_rel> = H.relativize(r)
from_rel, to_rel = Hrel.structure()
K = Hrel.base_field()

# elements in H for primality proving
def p_k(k):
    return 1 - r^k*a

# pk mod I
def p_k_mod(k, I):
    OI = H.quotient_ring(I)
    R = OI(r)
    return 1 - R^k * a

# odd part of N_{K|Q}(1 + r^k)
def c_k(k):
    ck = (K(1 + (r_rel)^k)).norm()
    while (ck % 2) == 0:
        ck //= 2
    return ZZ(ck)

# the j-invariant with CM by Z[r] over Q(a)
auto_Qa = [auto for auto in H.automorphisms() if auto(a) == a and (not auto(r) == r)][0]
for j in HT(hilbert_class_polynomial(-31)).roots(multiplicities=False):
    if auto_Qa(j) == j:
        jD = j
        break

# an elliptic curve of j-invariant jD 
E = EllipticCurve(j=jD)
A = E.a4()
B = E.a6()
# reduce factors of the coefficient
for l, e in factor(H.ideal(A)):
    while e >= 2:
        e -= 2
        if (l^3).divides(B):
            A //= l.gen(0)^2
            B //= l.gen(0)^3
# make B square in H
for l, e in factor(B):
    if not e % 2 == 0:
        A *= l^2
        B *= l^3
u = B // prod([l^e for l, e in factor(B)])
A *= u^2
B *= u^3
E = EllipticCurve(H, [A, B])
Bsqrt = sqrt(B)

assert B.is_square()
assert E.j_invariant() == jD
