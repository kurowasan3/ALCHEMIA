# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
ITEM = '/home/claude/chem2/out/items/THE-ELC-030/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

e = R('5.00')*R('1.93e4')/R('9.65e4'); print('電子', e); assert e == 1
O2 = R('0.560')/R('22.4'); e_el = 4*O2; e_ch = e - e_el
eff = e_ch/e*100; H2 = e_el/2; VH2 = H2*R('22.4')
print('O2', O2, 'eff', float(eff), 'H2 L', float(VH2)); assert eff == 90 and float(VH2) == 1.12
acid = R('0.500') + e_ch; water = 250 - (e_ch + H2)*R('18.0')
m = acid/(water/1000); print('acid', acid, 'water', float(water), 'molality', float(m))
assert round(float(m), 2) == 6.01
print('全検証OK:', ITEM)
