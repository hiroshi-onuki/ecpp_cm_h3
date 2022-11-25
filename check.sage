load("n31_def_symbols.sage")

for k in range(1, 10):
    Fk = ZZ(p_k(k).norm())
    Nk = (K(1 + (r_rel)^k)).norm()
    print(Fk)
    print(16*Nk^2)

