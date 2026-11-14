# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/THE-MOL-043/item.json'
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

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
assert check_equation('2CaCO3 + 2SO2 + O2 -> 2CaSO4 + 2CO2') is None
nC, nH2O, nS = 800/12, 50/1.0/2, 20/32
nO2 = nC + nH2O/2 + nS; assert sig(nO2, 2) == 80
assert sig(nO2*5*22.4/1000, 2) == 8.9
N2 = nO2*4; tot = nC + nH2O + nS + N2; dry = tot - nH2O
assert sig(nC/tot*100, 2) == 16 and sig(nC/dry*100, 2) == 17
ppm = nS/dry*1e6; print(ppm); assert sig(ppm, 2) == 1.6e3
r = nS*0.80; assert abs(r*100 - 50) < 1e-9 and abs(r*136 - 68) < 1e-9

print('ALL OK')
