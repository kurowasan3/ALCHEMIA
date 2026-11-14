ITEM = '/home/claude/chemk/outh/items/THE-THC-073/item.json'
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

import sympy as sp
T = sp.symbols('T')
Hf = {'C':0,'O2':0,'CO':-110.5,'CO2':-393.5,'Cr':0,'Cr2O3':-1139.7}
S = {'C':5.7,'O2':205.1,'CO':197.7,'CO2':213.8,'Cr':23.8,'Cr2O3':81.2}
def line(react, prod):
    H = sum(n*Hf[f] for f, n in prod) - sum(n*Hf[f] for f, n in react)
    s = sum(n*S[f] for f, n in prod) - sum(n*S[f] for f, n in react)
    return H, s/1000
a = line([('C',2),('O2',1)], [('CO',2)]); b = line([('C',1),('O2',1)], [('CO2',1)]); c = line([('Cr',sp.Rational(4,3)),('O2',1)], [('Cr2O3',sp.Rational(2,3))])
print(a, b, c)
assert abs(a[0]+221.0) < 1e-9 and abs(a[1]-0.1789) < 1e-9 and abs(b[1]-0.0030) < 1e-9 and abs(float(c[0])+759.8) < 1e-6 and abs(float(c[1])+0.1827) < 1e-6
G = lambda l: l[0] - T*l[1]
tab = float(sp.solve(sp.Eq(G(a), G(b)), T)[0]); tac = float(sp.solve(sp.Eq(G(a), G(c)), T)[0]); tbc = float(sp.solve(sp.Eq(G(b), G(c)), T)[0])
print(tab, tac, tbc); assert sig(tab,3) == 981 and sig(tac,3) == 1490 and sig(tbc,3) == 1970
# (a)-(c) が T>tac で負
assert float((G(a)-G(c)).subs(T, tac+1)) < 0 and float((G(a)-G(c)).subs(T, tac-1)) > 0
lo = -221.0 - (0.1789+0.180)*1200; hi = -221.0 - (0.1789+0.180)*1100; print(lo, hi)
assert sig(lo,3) == -652 and sig(hi,3) == -616
for e in ['2Cr2O3 + 6C -> 4Cr + 6CO', '2Cr2O3 + 3C -> 4Cr + 3CO2', '4Cr + 3O2 -> 2Cr2O3', 'C + CO2 -> 2CO']:
    assert check_equation(e) is None

print('全検証OK:', ITEM)
