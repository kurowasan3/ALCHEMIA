ITEM = '/home/claude/chemk/out/items/INO-ALK-026/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from fractions import Fraction as Fr
# 成分ごとの必要量（1e-3 mol/L）
w = {'Ca(HCO3)2': Fr(3,2), 'Mg(HCO3)2': Fr(1,2), 'MgCl2': Fr(4,5), 'CaSO4': Fr(1)}
lime = w['Ca(HCO3)2']*1 + w['Mg(HCO3)2']*2 + w['MgCl2']*1
soda = w['MgCl2']*1 + w['CaSO4']*1
print(lime, soda); assert lime == Fr(33,10) and soda == Fr(9,5)
gl = float(lime)*74; gs = float(soda)*106; print(gl, gs); assert sig(gl,3) == 244 and sig(gs,3) == 191
# Ca・Mg の収支：沈殿する Ca と加えた Ca＋元の Ca
Ca_in = w['Ca(HCO3)2'] + w['CaSO4'] + lime
CO3 = 2*w['Ca(HCO3)2'] + 2*w['Mg(HCO3)2'] + soda   # HCO3- から生じる CO3^2- と Na2CO3
assert Ca_in == CO3, (Ca_in, CO3)
OH = 2*lime; OH_used = 2*w['Ca(HCO3)2'] + 2*w['Mg(HCO3)2'] + 2*(w['Mg(HCO3)2'] + w['MgCl2'])
assert OH == OH_used
Na = 2*soda; Cl = 2*w['MgCl2']; SO4 = w['CaSO4']
assert Na == Cl + 2*SO4; print('Na+', float(Na)); assert sig(float(Na)*1e-3, 2) == 3.6e-3

print('全検証OK:', ITEM)
