# THE-THC-039 検証：分子とその陽イオンの結合エネルギー（熱化学サイクル）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['H2 -> 2H', 'H -> H+ + e-', 'H2 -> H2+ + e-', 'H2+ -> H + H+',
          'O2+ -> O + O+', 'N2+ -> N + N+', 'H2 -> H + H+ + e-']:
    r = check_equation(e); print(e, r or 'OK'); assert r is None
data = {'H': (436, 1312, 1488), 'N': (945, 1402, 1503), 'O': (498, 1314, 1165)}
res = {}
for k, (D, Ia, Im) in data.items():
    Dp = D + Ia - Im; res[k] = Dp
    print(k, 'D(M2+) =', Dp, 'stronger than M2?' , Dp > D)
assert res == {'H': 260, 'N': 844, 'O': 647}
# 文献値との比較（eV→kJ/mol, 96.485）：D0(H2+)=2.65 eV, D0(N2+)=8.71 eV, D0(O2+)=6.66 eV
for k, ev in (('H', 2.65), ('N', 8.71), ('O', 6.66)):
    print(k, 'lit', round(ev*96.485), 'calc', res[k])
# (3) 光子のエネルギー：E = 1.20e5 / λ(nm) kJ/mol
for lam in (85, 80, 68):
    E = 1.20e5/lam; print(lam, 'nm ->', round(E), 'kJ/mol', 'ionize H2' if E >= 1488 else '-', 'dissoc.ion.' if E >= 436+1312 else '-')
assert 1.20e5/85 < 1488 < 1.20e5/80 < 1748
lam_min = 1.20e5/1748; print('threshold dissociative ionization', round(lam_min, 1), 'nm')
assert round(lam_min, 1) == 68.6
# H2 の最低点と H2+ の最低点のエネルギー差 = I(H2)
print('OK')
