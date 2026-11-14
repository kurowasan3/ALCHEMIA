# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
from sympy import symbols, solve, Rational as Q
ITEM = '/home/claude/chem2/out/items/INO-NIT-015/item.json'
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

a, b, c, d, e = symbols('a b c d e')
sol = solve([3*a - d, 2*a - 4, 8*a + 2*b - 3*d - e, b - d, c - e], [a, b, c, d, e]); print(sol)
assert [sol[x] for x in (a, b, c, d, e)] == [2, 6, 10, 6, 10]
assert check_equation('2Ca3(PO4)2 + 6SiO2 + 10C -> 6CaSiO3 + 10CO + P4') is None
# 酸化数：PO4^3- 中の P
cc, q = parse_formula('PO4^3-'); oxP = q + 2*cc['O']; assert oxP == 5
cc, q = parse_formula('CO'); oxC = 2*cc['O']; assert oxC == 2
assert 4*oxP == 10*oxC == 20
n = Q(620, 310); P4 = n/2*124; C = n*5*12; CO = n*5*Q('22.4'); H3PO4 = n/2*4*98
print('P4 kg', P4, 'C kg', C, 'CO m3', CO, 'H3PO4 kg', H3PO4)
assert (P4, C, CO, H3PO4) == (124, 120, Q('224'), 392)
print('すべて OK')
