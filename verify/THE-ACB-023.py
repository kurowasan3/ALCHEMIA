# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from scipy.optimize import brentq
ITEM = '/home/claude/chem2/out/items/THE-ACB-023/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

Ka, Kb, Kw = 2.0e-5, 2.0e-5, 1.0e-14
# (1)
h = 10**-4.4; print('10^-4.4 =', h)
ratio = Ka / 4.0e-5; frac = ratio / (1 + ratio); V1 = 20.0 * frac
print('(1)', V1); assert round(V1, 1) == 6.7
# (2)
OH = sqrt(Kw / Ka * 0.050); pH2 = 14 + log10(OH); print('(2)', pH2); assert round(pH2, 1) == 8.7 and 8.0 < pH2 < 9.8
# (3)
r3 = (Kw / Kb) / 1e-8; f3 = 1 / (1 + r3); print('(3)', f3); assert round(f3 * 100) == 95
# (4)
f4 = 1e3 / (1 + 1e3); print('(4)', f4); assert round(f4 * 100, 1) == 99.9

# 厳密解（電荷のつり合い）で近似を確認
def pH_acid(vb):
    V = 20.0 + vb; C = 2.0 / V; Na = 0.10 * vb / V
    f = lambda lh: 10**lh + Na - C * Ka / (Ka + 10**lh) - Kw / 10**lh
    return -brentq(f, -14, 0)
def pH_base(va):
    V = 20.0 + va; C = 2.0 / V; Cl = 0.10 * va / V; Kn = Kw / Kb
    f = lambda lh: 10**lh + C * 10**lh / (10**lh + Kn) - Cl - Kw / 10**lh
    return -brentq(f, -14, 0)
print('厳密 pH(6.67mL)=', pH_acid(20.0/3), ' 中和点', pH_acid(20.0), ' NH3 95.2%', pH_base(20.0*f3), ' NH3中和点', pH_base(20.0))
assert abs(pH_acid(20.0/3) - 4.4) < 0.01 and abs(pH_acid(20.0) - 8.70) < 0.01 and abs(pH_base(20.0*f3) - 8.0) < 0.01
assert round(pH_base(20.0), 1) == 5.3
print('全検証OK:', ITEM)
