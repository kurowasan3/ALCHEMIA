# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/INO-TRN-041/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
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
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

from sympy import symbols, solve
for t, n in ((1.0, 2.0e-3), (4.0, 4.0e-3), (9.0, 6.0e-3)): assert abs(n*n/t - 4.0e-6) < 1e-12
assert abs((1.0e-2)**2/4.0e-6 - 25) < 1e-9
x = symbols('x', positive=True)
for oF, oM in ((5.70/72, 4*5.18/232), (0.0792, 0.0893)):
    s = float(solve(x*oF + x/9*oM - 8.0e-5, x)[0])
    assert sig(s*1e4,2) == 9.0 and sig(s/9*1e4,2) == 1.0, s
s = float(solve(x*5.70/72 + x/9*4*5.18/232 - 8.0e-5, x)[0])
fe = s*5.70/72 + 3*(s/9)*5.18/232
assert sig(fe*56/7.87*1e4,2) == 5.5
for e in ['2Fe + O2 -> 2FeO', '3Fe + 2O2 -> Fe3O4', '6FeO + O2 -> 2Fe3O4']: assert check_equation(e) is None

print('ALL OK')
