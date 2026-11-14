# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
ITEM = '/home/claude/chem2/outh/items/THE-RDX-023/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'\{a\+\}', b):
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

E = {'Fe2/Fe': R('-0.44'), 'Fe3/Fe2': R('0.77'), 'I2/I': R('0.54'), 'Sn2/Sn': R('-0.14'), 'Sn4/Sn2': R('0.15')}
EFe3Fe = (1*E['Fe3/Fe2'] + 2*E['Fe2/Fe'])/3; print('(1)', float(EFe3Fe)); assert round(float(EFe3Fe), 2) == -0.04
assert abs(float(EFe3Fe) - (-0.037)) < 0.001   # 表の値 −0.037 V
emf2 = E['Fe3/Fe2'] - E['Fe2/Fe']; print('(2)', float(emf2)); assert float(emf2) == 1.21
disp = E['Fe2/Fe'] - E['Fe3/Fe2']; assert disp < 0
e31 = E['Fe3/Fe2'] - E['I2/I']; e32 = E['I2/I'] - E['Fe3/Fe2']; eFeI = E['I2/I'] - E['Fe2/Fe']
print('(3)', float(e31), float(e32), float(eFeI)); assert e31 > 0 and e32 < 0 and float(eFeI) == 0.98
ESn = (2*E['Sn4/Sn2'] + 2*E['Sn2/Sn'])/4; print('(4)', float(ESn)); assert float(ESn) == 0.005
dSn = E['Sn2/Sn'] - E['Sn4/Sn2']; print('Sn2+ 不均化', float(dSn)); assert float(dSn) == -0.29
print('全検証OK:', ITEM)
