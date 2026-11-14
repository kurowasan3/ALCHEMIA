# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
from sympy import symbols, Rational as Q, simplify, solve
ITEM = '/home/claude/chem2/out/items/INO-CSI-010/item.json'
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

s = symbols('s')
nO = (4 - s) + s/2
charge = 4 - 2*nO
assert simplify(nO - (4 - s/2)) == 0 and simplify(charge + (4 - s)) == 0
for k in range(5): print('s =', k, ' O/Si =', nO.subs(s, k), ' 電荷/Si =', charge.subs(s, k))
ox = {'Si': 4, 'O': -2, 'H': 1, 'Ca': 2, 'Mg': 2, 'Be': 2, 'Al': 3, 'Na': 1}
for f, s_exp in [('CaMgSi2O6', 2), ('Be3Al2Si6O18', 2), ('Mg3Si4O10(OH)2', 3)]:
    c, q = parse_formula(f)
    tot = sum(ox[e]*n for e, n in c.items())
    oSil = c['O'] - c.get('H', 0)          # OH- の O を除く
    sval = solve(nO - Q(oSil, c['Si']), s)[0]
    print(f, dict(c), '酸化数の和', tot, ' s =', sval)
    assert tot == 0 and sval == s_exp
c, q = parse_formula('Si6O18^12-'); assert q == -12 and c['O'] == 6 * nO.subs(s, 2)
c, q = parse_formula('Na2SiO3'); assert solve(nO - Q(c['O'], 1), s)[0] == 2
print('すべて OK')
