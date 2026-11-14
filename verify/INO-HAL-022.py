# INO-HAL-022 verification: chlorine water with HCl impurity
from sympy import symbols, nsolve, solve, Rational as R, N, log, sqrt
import math
c = R(5, 100) * R(400, 100) / 20; cd = R(5, 100) * R(640, 100) / 20
d = cd - c
print('c', c, 'd', d); assert c == R(1, 100) and d == R(6, 1000)
H = 10 ** -1.885; print('H', H); assert round(H, 4) == 0.0130
H = 0.0130; x = H - float(d); K = H * H * x / (float(c) - x)
print('x', x, 'K', K); assert round(x, 4) == 0.0070 and round(K * 1e4, 1) == 3.9
print('ClO-/HClO', 3.0e-8 / H)
X = symbols('X', positive=True)
Kv = R(4, 10000)
r = [s for s in solve(X**3 - Kv * (c - X), X) if s.is_real and s > 0][0]
print('pure x', N(r), 'frac', N(r / c), 'pH', N(-log(r, 10)))
assert round(float(r / c), 2) == 0.85 and round(float(-log(r, 10)), 2) == 2.07
f = lambda v: v**3 - 4e-4 * (0.01 - v)
assert f(8.0e-3) < 0 < f(8.5e-3)
# first equivalence point
V = 26.40
Cl = 0.0160 * 20 / V; HClO = 0.0100 * 20 / V
h = math.sqrt(3.0e-8 * HClO)
cl2 = h * Cl * HClO / 4.0e-4
print('Cl-', Cl, 'HClO', HClO, 'H', h, 'pH', -math.log10(h), 'Cl2', cl2, cl2 / HClO)
assert round(-math.log10(h), 2) == 4.82 and round(cl2 * 1e6, 1) == 3.5
assert abs(math.log10(8.5) - 0.93) < 0.005 and abs(math.log10(1.51) - 0.18) < 0.005 and abs(math.sqrt(2.27) - 1.51) < 0.005
print('OK')
