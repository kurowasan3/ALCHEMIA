"""THE-ELC-042 検証：電気透析で脱塩される酢酸緩衝液の pH と近似の破れ（SymPy で厳密解と比較）"""
from sympy import Rational as R, symbols, solve, log, sqrt, N, nsolve
F = 96500; Ka = R(2, 10**5); Kw = R(1, 10**14)
ca0, cs0 = R(1, 10), R(1, 10)          # CH3COOH, CH3COONa 各 0.100 mol（1.00 L）
h = symbols('h', positive=True)
def exact_h(cs):
    # 電荷収支 [Na+]+[H+]=[A-]+[OH-]，物質収支 [HA]+[A-]=ca0+cs（酸 ca0 は膜を通らず一定）
    A = cs + h - Kw / h
    eq = h * A - Ka * (ca0 + cs - A)
    return nsolve(eq, h, 1e-4, prec=30)
def hh(cs): return Ka * ca0 / cs     # Henderson（[A-]=cs, [HA]=ca0）
# (1) 1.93 A × 2500 s
x1 = R(193, 100) * 2500 / F; print('電子 mol', x1); assert x1 == R(5, 100)
pH0 = -log(hh(cs0), 10); pH1 = -log(hh(cs0 - x1), 10)
print('pH 初め', float(N(pH0)), '後', float(N(pH1)), '厳密', float(-log(exact_h(cs0 - x1), 10)))
assert round(float(N(pH0)), 2) == 4.70 and round(float(N(pH1)), 2) == 4.40
# (2) pH 3.70 となる時間（Henderson）と厳密値
cs2 = Ka * ca0 / R(2, 10**4); x2 = cs0 - cs2; t2 = x2 * F / R(193, 100)
print('cs', cs2, 'x', x2, 't', float(t2)); assert cs2 == R(1, 100) and t2 == 4500
hx = R(2, 10**4); cs2e = Ka * (ca0 - hx) / hx - hx   # 厳密：[A-]=cs+h，[HA]=ca0-h（OH- 無視）
t2e = (cs0 - cs2e) * F / R(193, 100); print('厳密 t', float(t2e)); assert abs(t2e - t2) / t2 < 0.005
print('(2) の状態で [H+]/[Na+] =', float(hx / cs2))
# (3) 電子 0.0995 mol：cs = 5.0e-4
cs3 = cs0 - R(995, 10000); print('cs3', cs3)
hH = hh(cs3); print('Henderson [H+]', float(hH), '> cs3 なので矛盾'); assert hH > cs3
hq = solve(h**2 + (cs3 + Ka) * h - Ka * ca0, h)[0]; print('二次方程式', float(hq)); assert round(float(hq), 5) == 0.00118
print('√8.27e-6 を 2.88e-3 とした値', (2.88e-3 - 5.2e-4) / 2)
he = exact_h(cs3); print('厳密（OH- 含む）', float(he), 'pH', float(-log(he, 10)))
assert abs(float(he) - float(hq)) < 1e-8 and round(float(-log(he, 10)), 1) == 2.9
hp = sqrt(Ka * ca0); print('純粋な 0.100 M 酢酸の [H+]（近似）', float(hp))
print('(3) の状態で [H+]/[Na+] =', float(hq / cs3))
assert hq / cs3 > 2
print('THE-ELC-042 all OK')
