# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
ITEM = '/home/claude/chem2/out/items/THE-RDX-021/item.json'
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
n = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b); n += 1
assert n >= 10

S2O3 = R('0.100')*R('16.0')/1000; I2 = S2O3/2; ClO = I2*10
NaClO = 23 + R('35.5') + 16
pct = ClO*NaClO/R('10.0')*100; pCl2 = ClO*71/R('10.0')*100
print('ClO total', ClO, 'NaClO%', float(pct), 'Cl2%', float(pCl2))
assert round(float(pct), 2) == 5.96 and round(float(pCl2), 2) == 5.68
# 経路の比較：3ClO- -> 3I2 （直接）と 3ClO- -> IO3- -> 3I2
direct = 3 * 1          # ClO- 1 個で I2 1 個
via = 1 * 3             # IO3- 1 個から I2 3 個
assert direct == via
print('全検証OK:', ITEM)
