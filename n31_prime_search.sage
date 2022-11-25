import json

load("n31_def_symbols.sage")
load("tools.sage")
load("primality_proof.sage")

with open("n31_data/T.json", "r") as fp:
    T = json.load(fp)

a_mod_Fk = lambda r_Fk, k: r_Fk^(-k)
minpoly_a = lambda t: t^3 + t + 1
D = -31

threshold = 100
for k in range(2, 20000):
    if (k % T["period"]) in T["values"]:
        Fk = ZZ(p_k(k).norm())
        result = primality_proof(E, D, k, Fk, a_mod_Fk, minpoly_a)
        if result:
            print("F_%d is prime." % k)
        if k >= threshold:
            print(k)
            threshold += 100