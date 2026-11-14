# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/INO-TRN-032/item.json'
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
assert check_equation('Cu2CO3(OH)2 + 4HNO3 -> 2Cu(NO3)2 + CO2 + 3H2O') is None
assert check_equation('Cu4SO4(OH)6 + 6HNO3 -> 3Cu(NO3)2 + CuSO4 + 6H2O') is None
assert check_equation('Ba^2+ + SO4^2- -> BaSO4') is None
H, C, O, S, Cu, Ba = 1.0, 12, 16, 32, 63.5, 137
Mm = 2*Cu + C + 5*O + 2*H; Mb = 4*Cu + S + 10*O + 6*H; MB = Ba + S + 4*O
assert (Mm, Mb, MB) == (221.0, 452.0, 233)
nm = 0.224/22.4; nb = 1.165/MB
assert abs(nm*Mm + nb*Mb - 4.47) < 1e-9
assert sig(nb*Mb/4.47*100, 3) == 50.6
assert sig((2*nm + 4*nb)*Cu/4.47*100, 3) == 56.8
import json
ions = open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8').read()
assert 'Cu(NH3)4' in ions or 'Cu(NH_3)_4' in ions or '[Cu(NH3)4]' in ions

print('ALL OK')
