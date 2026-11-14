# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as Q
ITEM = '/home/claude/chem2/out/items/INO-AMP-015/item.json'
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

zn, mg, al = symbols('zn mg al')
s = solve([2*zn + 2*mg - 30, 2*zn - 10], [zn, mg]); print(s)
assert s[zn] == 5 and s[mg] == 10
H, O, Mg, Zn = Q(1), Q(16), Q('24.3'), Q('65.4')
mmax = s[zn]*(Zn + 2*O + 2*H)/1000 + s[mg]*(Mg + 2*O + 2*H)/1000
mend = s[mg]*(Mg + 2*O + 2*H)/1000
print('max', float(mmax), 'end', float(mend))
assert round(float(mmax), 2) == 1.08 and round(float(mend), 3) == 0.583
# Al と仮定
s2 = solve([3*al + 2*mg - 30, al - 10], [al, mg]); print('Al 仮定', s2); assert s2[mg] == 0
# 沈殿の物質量の合計：0→30 mL で 0→15 mmol（直線），40 mL で 10 mmol（図の高さ 180px:120px と一致）
assert s[zn] + s[mg] == 15 and s[mg] == 10
# 図鑑との照合
R = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
assert R['Zn^2+']['NaOH-ex']['kind'] == 'dissolve' and R['Mg^2+']['NaOH-ex']['kind'] == 'remain'
assert R['Zn^2+']['NH3-ex']['kind'] == 'dissolve' and R['Mg^2+']['NH3-ex']['kind'] == 'remain'
assert R['Al^3+']['NaOH-ex']['kind'] == 'dissolve'
print('すべて OK')
