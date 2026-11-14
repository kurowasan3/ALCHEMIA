ITEM = '/home/claude/chemx/x1/out/items/THE-GAS-063/item.json'
# 検証スクリプト（x1 g01）  python3 -I このファイル
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
assert nq == na or (nq == 0 and na <= 1), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

from fractions import Fraction as Fr
for e in ['C4H8 + 6O2 -> 4CO2 + 4H2O','2C4H10 + 13O2 -> 8CO2 + 10H2O']:
    assert check_equation(e) is None
# n from ratio
for n in range(2,8):
    nA = Fr(22,100) + Fr(2,100)*(Fr(n,2)-1)
    if nA/Fr(2,10) == Fr(6,5): print('n=',n); assert n==4
P = 0.440*8.31e3*400/20.0; chk('P',P,7.31e4,3)
chk('pH2O',P/3,2.44e4,3); assert P/3 < 2.5e5
assert 0.020*6 <= 0.200 and 0.020*6.5 <= 0.200
nA = Fr(22,100)+Fr(2,100)*Fr(3,2); r = nA/Fr(2,10); print(r); assert r==Fr(5,4)

print('ALL OK')
