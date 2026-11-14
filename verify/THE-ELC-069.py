# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/THE-ELC-069/item.json'
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

F = 9.65e4
Q = 0.100 * (64*60 + 20); n = Q / F
print(Q, n); assert abs(n - 4.00e-3) < 1e-9
cands = {'Mg2+': 24/2, 'Fe2+': 56/2, 'Fe3+': 56/3, 'Ni2+': 59/2, 'Zn2+': 65/2, 'Sn2+': 119/2, 'Sn4+': 119/4, 'Pb2+': 207/2, 'Cu2+': 64/2, 'Ag+': 108}
mX = 0.130 / n; mY = 0.432 / n
print(mX, mY)
mx = [k for k, v in cands.items() if abs(v - mX) < 0.15]; my = [k for k, v in cands.items() if abs(v - mY) < 0.15]
print(mx, my); assert mx == ['Zn2+'] and my == ['Ag+']
ag = 0.100*0.100 - n; print(ag/0.100); assert round(ag/0.100, 5) == 0.06
assert round(n*32*1000) == 128
print('OK')
