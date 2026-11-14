# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, solve
ITEM = '/home/claude/chem2/outh/items/THE-ELC-035/item.json'
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

Cu, Ni = Q('63.5'), Q('58.7')
nCu_dep = Q(635000) / Cu                     # 1.00e4 mol
ncu, nni = Q('97.8') / Cu, Q('1.80') / Ni
f = nni / (ncu + nni)
Ni_day = nCu_dep * f; dCu = nCu_dep * (1 - f) - nCu_dep
print('Ni/day', float(Ni_day), 'dCu/day', float(dCu))
assert round(float(Ni_day)) == 195 and round(float(dCu)) == -195
V = 40000
days = Q('0.200') * V / Ni_day; print('days', float(days)); assert round(float(days), 1) == 41.0
v, c0 = symbols('v c0')
s = solve([Q('0.200')*v - Ni_day, c0*v - Q('0.700')*v - Ni_day], [v, c0], dict=True)[0]
print('v', float(s[v]), 'c0', float(s[c0])); assert round(float(s[v])) == 976 and s[c0] == Q('0.9')
# 日ごとのシミュレーション（抜き取りなし → 和が一定）
cu, ni = 0.700 * V, 0.0
for d in range(41):
    cu += float(dCu); ni += float(Ni_day)
    assert abs((cu + ni) / V - 0.700) < 1e-12
print('41日後 Cu', cu / V, 'Ni', ni / V)
# 抜き取り・補給あり（毎日 電解→抜き取り v→補給 v）
cu, ni = 0.700 * V, 0.0; vv = float(s[v]); cc = float(s[c0])
for d in range(2000):
    cu += float(dCu); ni += float(Ni_day)
    cu -= vv * cu / V; ni -= vv * ni / V
    cu += vv * cc
print('定常', cu / V, ni / V)
# 収束先（電解直後に抜き取る順序なので抜き取り時点の濃度が 0.700/0.200）
cu_b, ni_b = (cu + float(dCu)) / V, (ni + float(Ni_day)) / V
print('抜き取り時点', cu_b, ni_b); assert abs(cu_b - 0.700) < 1e-6 and abs(ni_b - 0.200) < 1e-6
# Pb を含む場合：Pb x mol → e 2x，SO4 −x，陰極 Cu2+ −x（陽極から Cu2+ 補給なし）→ 2価陽イオン −x, 電荷 2(−x) − 2(−x) = 0
x = 1; assert 2*(-x) - 2*(-x) == 0
print('すべて OK')
