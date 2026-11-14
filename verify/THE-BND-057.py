# THE-BND-057 検証：気体の原子対で電子が移る距離（イオン対と中性原子のエネルギーの交差）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['K + Br -> K+ + Br-', 'KBr -> K + Br', 'KBr -> K+ + Br-']:
    r = check_equation(e); print(e, r or 'OK'); assert r is None
# k = N_A e^2/(4πε0) in kJ·nm/mol
NA, e, eps0 = 6.02214076e23, 1.602176634e-19, 8.8541878128e-12
import math
k = NA*e**2/(4*math.pi*eps0)*1e9/1000
print('k =', k, 'kJ nm/mol'); assert abs(k - 138.9) < 0.1
k = 1.39e2
I = {'Li': 520, 'Na': 496, 'K': 419, 'Cs': 376}
EA = {'F': 328, 'Cl': 349, 'Br': 325, 'I': 295}
pairs = [('Cs', 'Cl'), ('K', 'Br'), ('Li', 'F'), ('Na', 'I')]
rc = {}
for m, x in pairs:
    d = I[m] - EA[x]; rc[m+x] = k/d
    print(m+x, 'I-EA =', d, 'r_c =', round(k/d, 3), 'nm')
assert round(rc['KBr'], 2) == 1.48
order = sorted(rc, key=rc.get, reverse=True); print(order)
assert order == ['CsCl', 'KBr', 'LiF', 'NaI']
# (3) KBr 分子 r = 0.282 nm
E = -k/0.282 + (I['K'] - EA['Br']); print('E(model) =', round(E, 1), 'kJ/mol ->', round(-E))
assert round(-E, -1) == 400
# (4) 無限遠でイオンの方が安定な組合せはあるか：min I > max EA
print('min I', min(I.values()), 'max EA', max(EA.values())); assert min(I.values()) > max(EA.values())
print('OK')
