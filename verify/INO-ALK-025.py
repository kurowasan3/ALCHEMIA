ITEM = '/home/claude/chemk/out/items/INO-ALK-025/item.json'
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

from collections import Counter
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import parse_formula
M = {'H':1.0,'C':12,'N':14,'O':16,'Na':23,'Cl':35.5,'Ca':40}
fm = lambda f: sum(M[e]*n for e, n in parse_formula(f)[0].items())
def side(s):
    c = Counter()
    for t in s.split(' + '):
        t = t.strip(); n = 1
        mm = re.match(r'^(\d+)(.*)$', t)
        if mm: n, t = int(mm.group(1)), mm.group(2)
        c[t] += n
    return c
def net(eqs):
    tot = Counter()
    for coef, e in eqs:
        L, R = e.split(' -> ')
        for f, n in side(L).items(): tot[f] -= coef*n
        for f, n in side(R).items(): tot[f] += coef*n
    return {f: n for f, n in tot.items() if n}
a='NaCl + NH3 + CO2 + H2O -> NaHCO3 + NH4Cl'; b='2NaHCO3 -> Na2CO3 + CO2 + H2O'; c='CaCO3 -> CaO + CO2'; d='CaO + H2O -> Ca(OH)2'; e='2NH4Cl + Ca(OH)2 -> CaCl2 + 2NH3 + 2H2O'
n1 = net([(2,a),(1,b),(1,c),(1,d),(1,e)]); print(n1)
assert n1 == {'NaCl':-2,'CaCO3':-1,'Na2CO3':1,'CaCl2':1}
n3 = net([(2,a),(1,b)]); print(n3)
assert n3 == {'NaCl':-2,'NH3':-2,'CO2':-1,'H2O':-1,'Na2CO3':1,'NH4Cl':2}
ae1 = fm('Na2CO3')/(fm('Na2CO3')+fm('CaCl2'))*100; print(ae1); assert sig(ae1,3) == 48.8
ae2 = fm('Na2CO3')/(fm('Na2CO3')+2*fm('NH4Cl'))*100; print(ae2); assert sig(ae2,3) == 49.8
n = 1e6/fm('Na2CO3')
print(n*fm('CaCl2')/1e6, 2*n*fm('NH4Cl')/1e6, 2*n*fm('NaCl')/1e6)
assert sig(n*fm('CaCl2')/1e6,3) == 1.05 and sig(2*n*fm('NH4Cl')/1e6,3) == 1.01 and sig(2*n*fm('NaCl')/1e6,3) == 1.10

print('全検証OK:', ITEM)
