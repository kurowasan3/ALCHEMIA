# INO-PPT-031 検証：PbCl2 の溶解度積と第1属の分離の不完全さ
import json, re, sys
from sympy import Rational as Q, N, sqrt, cbrt
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
Ksp25, Ksp100, KspAg, KspPbS, S2 = Q(16, 10**6), Q(4, 10**3), Q(18, 10**11), Q(3, 10**28), Q(1, 10**22)
V = Q(1, 10); Pb0, Ag0, Cl = Q(2, 100), Q(1, 100), Q(1, 10)
Pb1 = Ksp25 / Cl**2; Ag1 = KspAg / Cl
print('[Pb2+] left', N(Pb1, 3), 'fraction', N(Pb1 / Pb0, 3), '[Ag+] left', N(Ag1, 3))
assert Pb1 == Q(16, 10**4) and Pb1 / Pb0 == Q(8, 100)
print('PbS Q', N(Pb1 * S2, 3), '>', N(KspPbS)); assert Pb1 * S2 > KspPbS
# (2) 0.1 % 以下にする [Cl-]
Cl_need = sqrt(Ksp25 / (Pb0 / 1000)); print('Cl needed', N(Cl_need, 3)); assert round(float(Cl_need), 2) == 0.89
# (3) 希薄な Pb2+
Pb_d = Q(1, 1000); print('Q dilute', N(Pb_d * Cl**2), '< Ksp'); assert Pb_d * Cl**2 < Ksp25
Pb_min = Ksp25 / Cl**2; assert Pb_min == Q(16, 10**4)
# (4) 熱湯
nppt = (Pb0 - Pb1) * V; s100 = cbrt(Ksp100 / 4); s25 = cbrt(Ksp25 / 4)
Vhot = nppt / s100; print('n ppt', N(nppt, 3), 's100', N(s100, 3), 'V hot mL', N(Vhot * 1000, 3))
assert nppt == Q(184, 10**5) and s100 == Q(1, 10) and round(float(Vhot * 1000)) == 18
remain = s25 * Vhot; frac = 1 - remain / nppt
print('s25', N(s25, 4), 'crystallised fraction', N(frac, 3)); assert round(float(s25), 4) == 0.0159 and round(float(frac) * 100) == 84
it = json.load(open('/home/claude/chem/out2/items/INO-PPT-031/item.json')); k = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') + eq.count('<=>') == 1 and 'boxed' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); k += 1; print('eq ok', eq)
assert k >= 2
print('OK')
