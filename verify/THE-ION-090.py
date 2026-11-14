# THE-ION-090 検証：PbSO4 の沈殿で Pb2+ を 1.0e-5 mol/L 以下にする Na2SO4 水溶液の最小体積
from sympy import symbols, solve, Rational, nsimplify
Ksp = 1.6e-8; nPb = 0.020*0.100; target = 1.0e-5; V0 = 0.100
def vmin(c):
    V = symbols('V', positive=True)
    S = Ksp/target                     # 必要な [SO4^2-]
    # 加えた SO4 - 沈殿した SO4（= nPb - 残る Pb） = 残る SO4
    sol = solve(c*V - (nPb - target*(V0+V)) - S*(V0+V), V)
    return [float(s) for s in sol]
for c in (0.10, 0.010, 0.0015):
    print(c, vmin(c))
print('必要[SO4]=', Ksp/target)
v1 = vmin(0.10)[0]*1000; v2 = vmin(0.010)[0]*1000
print(f'(1) {v1:.2f} mL -> 22 mL  (2) {v2:.1f} mL -> 2.6e2 mL')
assert round(v1) == 22 and round(v2, -1) == 260
assert vmin(0.0015) == []   # 解なし（いくら加えても届かない）
# 0.0015 mol/L で無限に加えた極限の [SO4] は 0.0015 < 1.6e-3
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['Pb^2+ + SO4^2- -> PbSO4', 'PbSO4 <=> Pb^2+ + SO4^2-']:
    assert check_equation(e) is None, e
print('OK')
