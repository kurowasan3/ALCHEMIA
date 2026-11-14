# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/INO-ALK-022/item.json'
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

from sympy import symbols, solve, Rational
x, y = symbols('x y')
s = solve([100*x + 197*y - Rational(1985,1000), 44*x - Rational(440,1000)], [x, y])
assert s[x] == Rational(1,100) and s[y] == Rational(1,200)
loss = 0.440*0.90; xc = loss/44; yb = (1.985 - xc*100)/197
assert sig(xc/0.1,3) == 0.0900 and sig(yb/0.1,3) == 0.0551, yb
# (4) Ba の一部分解
loss2 = 0.440 + 0.1*0.005*44; xc2 = loss2/44; yb2 = (1.985 - xc2*100)/197; assert xc2 > 0.01 and yb2 < 0.005
assert sig(100/197*100,3) == 50.8
assert check_equation('BaCO3 -> BaO + CO2') is None

print('ALL OK')
