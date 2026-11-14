# THE-BND-098 検証：分子振動の数・換算質量による同位体置換の振動数比
from sympy import Rational as R, sqrt, N
mods = {'CS2': (3, True), 'O3': (3, False), 'HCN': (3, True), 'NH3': (4, False), 'SF6': (7, False)}
for k, (n, lin) in mods.items():
    print(k, 3*n - (5 if lin else 6))
assert [3*n - (5 if lin else 6) for n, lin in mods.values()] == [4, 3, 4, 6, 15]
def mu(a, b): return a*b/(a + b)
rOH = sqrt(mu(16, 2)/mu(16, 1)); rSH = sqrt(mu(32, 2)/mu(32, 1)); rCH = sqrt(mu(12, 2)/mu(12, 1))
print('O-H/O-D', N(rOH), 'S-H/S-D', N(rSH), 'C-H/C-D', N(rCH), 'limit', N(sqrt(2)))
assert f'{float(rOH):.3g}' == '1.37' and f'{float(rSH):.3g}' == '1.39'
od = 3660/float(rOH)
print('O-D wavenumber', od)
assert f'{od:.2g}' == '2.7e+03'
print(round(od, -1))
