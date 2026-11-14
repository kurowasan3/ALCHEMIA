# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/ORG-DET-028/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
print('小問の数', np_)

from math import log10, floor, sqrt
def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import mpmath
Ka = 3.0e-10; K1 = 4.5e-7; K2 = 4.7e-11; Kw = 1.0e-14; S0 = 5.0e-3
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
m = Chem.MolFromSmiles('Oc1ccc2ccccc2c1'); assert CalcMolFormula(m) == 'C10H8O'
M = 10*12 + 8*1.0 + 16; assert M == 144
K = Ka / K2; Kp = Ka / K1
print('(1)', K, Kp); assert sig(K, 2) == 6.4 and sig(Kp, 2) == 6.7e-4
# (2) 近似
C = 0.100
y = (-K*S0 + sqrt((K*S0)**2 + 4*K*S0*C)) / 2
tot2 = S0 + y; print('(2) 近似 y', y, '全濃度', tot2); assert sig(tot2, 2) == 0.048
# 厳密：電荷のつり合い Na+ + H+ = 2CO3 + HCO3 + OH + ArO（炭素の総量一定）
def full(Na, Ct):
    def f(L):
        h = 10**(-L)
        D = h*h + K1*h + K1*K2
        return Na + h - (Ct*(K1*h + 2*K1*K2)/D + Kw/h + Ka*S0/h)
    L = mpmath.findroot(f, 9.5)
    h = 10**(-L); D = h*h + K1*h + K1*K2
    return dict(pH=float(L), ArO=float(Ka*S0/h), OH=float(Kw/h), CO2=float(C*h*h/D), CO3=float(Ct*K1*K2/D), HCO3=float(Ct*K1*h/D))
r2 = full(0.200, C); print('(2) 厳密', r2)
assert sig(S0 + r2['ArO'], 2) == 0.048
assert r2['OH'] < 0.01*r2['ArO'] and r2['CO2'] < 0.01*r2['ArO']
assert round(r2['pH'], 1) == 10.5 or round(r2['pH'], 1) == 10.4
h_ap = K2*y/(C - y); print('近似の pH', -log10(h_ap)); assert round(-log10(h_ap), 1) == 10.5 or round(-log10(h_ap),1)==10.4
print('加える前の OH-', sqrt(C*Kw/K2))
# (3)
h3 = sqrt(K1*K2); A3 = Ka*S0/h3; tot3 = S0 + A3
r3 = full(0.100, C); print('(3) 近似', -log10(h3), A3, tot3, ' 厳密', r3, S0 + r3['ArO'])
assert sig(tot3, 2) == 5.3e-3 and sig(S0 + r3['ArO'], 2) == 5.3e-3
# (4)
n = 0.50 / M; cap2 = tot2*0.100; cap3 = tot3*0.100
print('(4) n', n, 'cap Na2CO3', cap2, 'cap NaHCO3', cap3, '残り g', (n - cap3)*M)
assert cap2 > n and cap3 < n
ppt = (n - S0*0.100) * M; print('析出 g', ppt); assert sig(ppt, 2) == 0.43
assert Ka/1e-8 <= 0.03 + 1e-12
Kco2 = K1/Ka; print('ArO- + CO2 の K', Kco2); assert sig(Kco2, 2) == 1.5e3
print('ORG-DET-028 OK')
