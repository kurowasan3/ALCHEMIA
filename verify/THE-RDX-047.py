ITEM = '/home/claude/chemx/tokodai/out/items/THE-RDX-047/item.json'
# 検証スクリプト（東工大 g04）  python3 -I この.py
import sys, json, re, math
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
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

from fractions import Fraction as F
c = F(2,5)*F(5,100)*25/25
chk('(1)', float(c), 2.00e-2, 3)
rest = F(5,2)*c*30
assert rest == F(3,2)
mno2 = F(5,100)*50 - rest
assert mno2 == 1
chk('(2)', float(mno2*87/1000/F(12,100)*100), 72.5, 3)
chk('(4)', float(F(12,100)/87*1000/F(5,100)), 27.6, 3)
for e in ['2KMnO4 + 5H2C2O4 + 3H2SO4 -> 2MnSO4 + K2SO4 + 10CO2 + 8H2O','MnO2 + 4HCl -> MnCl2 + Cl2 + 2H2O','MnO2 + H2C2O4 + H2SO4 -> MnSO4 + 2CO2 + 2H2O']:
    assert check_equation(e) is None

print('ALL OK')
