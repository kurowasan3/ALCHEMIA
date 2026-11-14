# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
ITEM = '/home/claude/chem2/out/items/THE-RDX-019/item.json'
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

cM, cO = R('2.00e-3'), R('5.00e-3')
eOx = 2*cO*R('10.0e-3')
eS = 5*cM*R('14.20e-3') - eOx
eB = 5*cM*R('10.20e-3') - eOx
dE = eS - eB
O2mg = dE/4*32*1000
COD = O2mg*10
print('(1)', float(eS), float(eB), float(dE), float(O2mg), float(COD)); assert round(float(COD), 1) == 3.2
cap = 5*cM*R('10.0e-3')/2
CODmax = cap/4*32*1000*10
print('(4)', float(cap), float(CODmax)); assert round(float(CODmax), 1) == 4.0
# 図鑑：AgCl
import json as J
ions = J.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))
s = J.dumps(ions, ensure_ascii=False)
assert 'AgCl' in s; print('図鑑に AgCl あり')
print('全検証OK:', ITEM)
