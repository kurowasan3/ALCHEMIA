# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from scipy.optimize import brentq
ITEM = '/home/claude/chem2/out/items/THE-ION-044/item.json'
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

Ka, Kw = 2.0e-5, 1.0e-14
pKa = -log10(Ka)
def exact(HA, A, Cl):   # 1 L, Na+ = A, 加えた HCl = Cl
    CT = HA + A
    f = lambda lh: 10**lh + A - Cl - CT*Ka/(Ka + 10**lh) - Kw/10**lh
    return -brentq(f, -14, 0)
sol = {'A': (0.10, 0.10), 'B': (0.050, 0.050), 'C': (0.16, 0.040), 'D': (0.020, 0.020)}
L2, L3, L7 = 0.30, 0.48, 0.85
given = {'A0': 4.70, 'C0': 4.70 - 2*L2, 'A': 4.70 + L2 - L3, 'B': 4.70 + L3 - L7, 'C': 4.70 - 2*L3, 'D': 0.5*(7 - 3*L2)}
print('与えた対数値', given)
assert round(given['A0'], 2) == 4.70 and round(given['C0'], 2) == 4.10 and round(given['A'], 2) == 4.52
assert round(given['B'], 2) == 4.33 and round(given['C'], 2) == 3.74 and round(given['D'], 2) == 3.05
for k, (ha, a) in sol.items():
    p0 = pKa + log10(a/ha)
    if k != 'D':
        p1 = pKa + log10((a - 0.020)/(ha + 0.020))
    else:
        p1 = -log10(sqrt(Ka*0.040))
    e0, e1 = exact(ha, a, 0.0), exact(ha, a, 0.020)
    print(k, '近似', round(p0, 3), round(p1, 3), '厳密', round(e0, 3), round(e1, 3), '変化', round(p0 - p1, 3))
    assert abs(p0 - e0) < 0.01 and abs(p1 - e1) < 0.01
    key = k if k != 'D' else 'D'
    assert abs(p1 - given[key]) < 0.015
# 変化の順序 A < C <= B
dA = 4.70 - given['A']; dB = 4.70 - given['B']; dC = given['C0'] - given['C']
print('ΔpH', dA, dB, dC); assert dA < dC <= dB + 0.02
print('全検証OK:', ITEM)
