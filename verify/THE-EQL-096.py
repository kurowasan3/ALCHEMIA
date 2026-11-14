"""THE-EQL-096 検証"""
from sympy import symbols, expand, solve, simplify, Rational as Q
k1, k2, R0, L0, RL, K = symbols('k1 k2 R0 L0 RL K', positive=True)
v = k1 * L0 * (R0 - RL) - k2 * RL
assert expand(v - (-(k1 * L0 + k2) * RL + k1 * L0 * R0)) == 0
RLeq = solve(v, RL)[0]
X = simplify((RLeq / R0).subs(k1, K * k2))
assert simplify(X - K * L0 / (1 + K * L0)) == 0
Ls = [1.0e-8, 2.0e-8, 3.0e-8]; al = [2.5e-3, 4.0e-3, 5.5e-3]
slope = (al[2] - al[0]) / (Ls[2] - Ls[0]); icpt = al[0] - slope * Ls[0]
assert all(abs(slope * L + icpt - a) < 1e-15 for L, a in zip(Ls, al))
Kv = slope / icpt
print(slope, icpt, Kv, 1 / Kv, 9 / Kv)
assert f'{slope:.2g}' == '1.5e+05' and f'{icpt:.2g}' == '0.001'
assert f'{1/Kv:.2g}' == '6.7e-09' and f'{9/Kv:.2g}' == '6e-08'
print('THE-EQL-096 OK')
