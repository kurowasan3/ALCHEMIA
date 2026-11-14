# THE-EQL-098 検証：2 段階で働く酵素（E と E* の往復）の定常状態
from sympy import symbols, solve, simplify, Rational as Q, Eq
k1, k2, S, B, ET, E, Es = symbols('k1 k2 S B E_T E Es', positive=True)
# 定常状態：k1[S][E] = k2[B][E*],  [E] + [E*] = E_T
sol = solve([Eq(k1*S*E, k2*B*Es), Eq(E + Es, ET)], [E, Es], dict=True)[0]
v = simplify(k1*S*sol[E]); print('(1) v =', v)
assert simplify(v - k1*k2*S*B*ET/(k1*S + k2*B)) == 0
print('    [E]/[E*] =', simplify(sol[E]/sol[Es]))
# (2) [S] が半分で v が 2.0 % 低下 → r = [E]/[E*] = k2[B]/(k1[S]) を求める
r = symbols('r', positive=True)
rs = solve(Eq((1 + r)/(1 + 2*r), Q(98, 100)), r)[0]; print('(2) r =', rs, float(rs))
assert rs == Q(1, 48) and round(float(rs), 3) == 0.021
print('    fraction of E* =', float(1/(1 + rs)))
# (3) v が最大値 k2[B]E_T の 90 % 以上でいられる [S] の範囲
Sfrac = rs/Q(1, 9); print('(3) [S]/[S]0 >=', Sfrac, float(Sfrac)); assert Sfrac == Q(3, 16)
print('    consumed fraction at that point =', float(1 - Sfrac))
# (4) [B] を 2 倍 vs [S] を 2 倍
fB = 2*(1 + rs)/(1 + 2*rs); fS = (1 + rs)/(1 + rs/2)
print('(4) B x2 ->', float(fB), ' S x2 ->', float(fS))
assert round(float(fB), 2) == 1.96 and round(float(fS), 2) == 1.01
print('OK')
