# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N
ITEM = '/home/claude/chem2/out/items/INO-TRN-023/item.json'
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

C, O, Si, Ca, Fe = 12, 16, 28, 40, 56
ore = 1000
n_SiO2 = Q(ore) * Q(120, 1000) / (Si + 2*O)
m_CaCO3 = n_SiO2 * (Ca + C + 3*O); m_slag = n_SiO2 * (Ca + Si + 3*O)
n_Fe2O3 = Q(ore) * Q(840, 1000) / (2*Fe + 3*O)
m_Fe = 2 * n_Fe2O3 * Fe
pig = m_Fe / Q(955, 1000)
m_C = 3 * n_Fe2O3 * C + pig * Q(45, 1000)
print('CaCO3', N(m_CaCO3), 'slag', N(m_slag), 'Fe', N(m_Fe), 'pig', N(pig, 6), 'C', N(m_C, 6))
assert m_CaCO3 == 200 and m_slag == 232 and m_Fe == 588
assert round(float(pig)) == 616 and round(float(m_C)) == 217
print('素朴な誤り 588*1.045 =', 588 * 1.045)
print('すべて OK')
