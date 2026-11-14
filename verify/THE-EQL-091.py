# THE-EQL-091 検証：触媒表面の2段階反応と中間体の定常状態
import sympy as sp
k1, k2 = sp.Float('3.0e6'), sp.Float('5.0e11')
A, B = sp.Float('4.0e-10'), sp.Float('2.0e-9')
AB = sp.symbols('AB', positive=True)
sol = sp.solve(sp.Eq(k1*A*B, k2*AB*B), AB)[0]
v = k2*sol*B
print('[AB]=', sol, ' v=', v)
assert f"{float(sol):.1e}" == '2.4e-15' and f"{float(v):.1e}" == '2.4e-12'
# [B] を2倍
sol2 = sp.solve(sp.Eq(k1*A*2*B, k2*AB*2*B), AB)[0]; v2 = k2*sol2*2*B
print('B x2:', sol2/sol, v2/v); assert sol2 == sol and abs(v2/v - 2) < 1e-12
# k2 を 1/100
sol3 = sp.solve(sp.Eq(k1*A*B, k2/100*AB*B), AB)[0]; v3 = k2/100*sol3*B
print('k2/100:', sol3/sol, v3/v); assert abs(sol3/sol - 100) < 1e-9 and abs(v3/v - 1) < 1e-12
# 1 cm^2, 1 s あたりの分子数
print('molecules', float(v*6.0e23))
