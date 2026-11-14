# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
from sympy import Rational as R
ITEM = '/home/claude/chem2/out/items/THE-RDX-020/item.json'
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

M = {'H': 1.0, 'O': 16, 'S': 32, 'Fe': 56}
mFeSO4 = 56 + 32 + 64 + 7*18.0; mFe2O3 = 2*56 + 3*16
print('式量', mFeSO4, mFe2O3); assert mFeSO4 == 278 and mFe2O3 == 160
Fe2 = 5*R('0.0200')*R('12.0')/1000 / R('0.0100'); print('(1)', Fe2); assert Fe2 == R('0.120')
FeT = R('0.240')/160*2 / R('0.0200'); Fe3 = FeT - Fe2; print('(2)', FeT, Fe3, Fe3/FeT)
assert FeT == R('0.150') and Fe3 == R('0.030') and Fe3/FeT == R(1, 5)
w = FeT*R('0.100')*278; O2 = Fe3*R('0.100')/4*R('22.4')*1000; print('(3)', float(w), float(O2))
assert round(float(w), 2) == 4.17 and round(float(O2), 1) == 16.8
ions = open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8').read()
for k in ('Fe(OH)2', 'Fe(OH)3'):
    i = ions.find('"f": "%s"' % k); print(k, ions[i:i + 60].replace('\n', ' '))
    assert i > 0
print('全検証OK:', ITEM)
