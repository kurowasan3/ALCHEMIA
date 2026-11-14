# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
from sympy import Rational as Q
from collections import Counter
ITEM = '/home/claude/chem2/out/items/INO-NIT-014/item.json'
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

# 酸化数
ox = {}
for f in ['NO3^-', 'NO2^-', 'NO', 'N2O', 'N2']:
    c, q = parse_formula(f)
    ox[f] = Q(q + 2*c.get('O', 0), c['N'])
print(ox); assert [ox[k] for k in ox] == [5, 3, 2, 1, 0]

# 全体式を半反応式から再構成（左辺 − 右辺 の化学種の係数）
def side(t):
    d = Counter()
    for term in t.split(' + '):
        m = re.match(r'(\d*)(.+)', term.strip()); d[m.group(2)] += int(m.group(1) or 1)
    return d
def eq(e):
    l, r = e.split(' -> '); d = Counter(); d.update(side(l)); d.subtract(side(r)); return d
red = eq('2NO3- + 12H+ + 10e- -> N2 + 6H2O'); oxi = eq('CH3OH + H2O -> CO2 + 6H+ + 6e-')
tot = Counter()
for k, v in red.items(): tot[k] += 3*v
for k, v in oxi.items(): tot[k] += 5*v
tot = {k: v for k, v in tot.items() if v}
print(tot); assert tot == dict(eq('6NO3- + 5CH3OH + 6H+ -> 3N2 + 5CO2 + 13H2O'))

nNO3 = Q('0.62') * 1000 / 62; assert nNO3 == 10
m4 = nNO3 * 5 / 6 * 32; m5 = (nNO3 * Q(8, 10) * 5 + nNO3 * Q(2, 10) * 4) / 6 * 32
print(float(m4), float(m5)); assert round(float(m4), -1) == 270 and round(float(m5), -1) == 260
print('すべて OK')
