import json

load("n31_def_symbols.sage")
load("tools.sage")
load("condition_k_mod.sage")

D = -31
cs = [[2, 3],[2, 5]] #norm(2r + 3) = 47, norm(2r + 5) = 67

T0 = calcT0(E)
with open("n31_data/T0.json", "w") as fp:
    json.dump(T0, fp)

T1 = calcT1(E, cs, D)
with open("n31_data/T1.json", "w") as fp:
    json.dump(T1, fp)

T2 = calcT2(E, r, D)
with open("n31_data/T2.json", "w") as fp:
    json.dump(T2, fp)

T = integrate_results([T0, T1, T2])
with open("n31_data/T.json", "w") as fp:
    json.dump(T, fp)
