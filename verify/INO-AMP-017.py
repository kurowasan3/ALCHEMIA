# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
ITEM = '/home/claude/chem2/out/items/INO-AMP-017/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; d = 1
        while d:
            if s[j] == '{': d += 1
            elif s[j] == '}': d -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = {b for k in ('problem', 'answer', 'solution') for b in ce_blocks(it[k]) if len(re.findall(r'->|<=>', b)) == 1}
for e in sorted(eqs):
    r = check_equation(e); assert r is None, (e, r); print('式OK:', e)

O, Cr, Pb = 16, 52, 207
PbO2, Pb3O4, PbO, PbCrO4 = Pb + 2*O, 3*Pb + 4*O, Pb + O, Pb + Cr + 4*O
assert (PbO2, Pb3O4, PbO, PbCrO4) == (239, 685, 223, 323)
n = Q('2.39') / PbO2
m34 = n * Pb3O4
pct = m34 / 10 * 100
VCl2 = n * Q('22.4')
nPb2 = 2*n + (10 - m34) / PbO
mppt = nPb2 * PbCrO4
print('n', float(n), 'Pb3O4 %', float(pct), 'Cl2 L', float(VCl2), 'Pb2+ mol', float(nPb2), 'PbCrO4 g', float(mppt))
assert n == Q(1, 100) and pct == Q('68.5') and VCl2 == Q('0.224') and round(float(mppt), 1) == 11.0
R = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
assert R['Pb^2+']['CrO4']['f'] == 'PbCrO4' and R['Pb^2+']['CrO4']['color'] == '黄'
print('すべて OK')
