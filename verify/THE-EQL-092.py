# THE-EQL-092 検証：N2O4 ⇄ 2NO2 の光路中の NO2 量（全圧と光路長の比較）
import sympy as sp
Kp = sp.Rational(15, 100)
x = sp.symbols('x', positive=True)
def pno2(P):
    s = [r for r in sp.solve(sp.Eq(x**2/(P - x), Kp), x) if 0 < r < P]
    assert len(s) == 1; return s[0]
pC, pD = pno2(sp.Rational(3, 100)), pno2(sp.Rational(1, 100))
print('pNO2 C, D', float(pC), float(pD))
assert f"{float(pC):.2g}" == '0.026' and f"{float(pD):.2g}" == '0.0094'
nC, nD = pC*15, pD*45
print('column C, D', float(nC), float(nD), 'ratio D/C', float(nD/nC)); assert round(float(nD/nC), 2) == 1.10
# 全分子数（P×L）は等しい
assert sp.Rational(3, 100)*15 == sp.Rational(1, 100)*45
# Ar を加えても体積一定なら分圧は不変 → Kp の式は同じ → 変化なし
print('OK')
