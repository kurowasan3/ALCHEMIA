# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/THE-SOL-095/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

import sympy as sp
R = 8.31e3; T = 300; V = 0.500; g = V/(R*T)
kC, kO, kD, pw = 1.40e-8, 1.30e-8, 3.30e-7, 4.00e3
gC, gO, gD = g+kC, g+kO, g+kD
assert round(gC*1e7, 2) == 2.15
def P2(a, b, gd=gD):
    nC = a*gC; nO = b*gO
    return (nO - 2*nC)/gO + nC/gd + pw
p2 = P2(2.00e4, 7.00e4); print('P2', p2); assert float('%.3g' % p2) == 4.19e4
a, b = sp.symbols('a b')
s = sp.solve([a + b + pw - 9.40e4, P2(a, b) - 4.19e4], [a, b]); print(s)
assert float('%.2g' % s[a]) == 2.0e4 and float('%.2g' % s[b]) == 7.0e4
s2 = sp.solve([a + b + pw - 9.40e4, P2(a, b, g) - 4.19e4], [a, b]); print('ignore CO2', s2)
assert float('%.2g' % s2[a]) == 2.7e4
assert round(kD/gD, 2) == 0.62
r = 2*gC/gO; print(r); assert round(r, 2) == 2.01
assert (7.00e4*gO - 2*2.00e4*gC) > 0
print('OK')
