ITEM = '/home/claude/chemx/x1/out/items/THE-ELC-074/item.json'
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

for e in ['2H2O -> O2 + 4H+ + 4e-','Cu^2+ + 2e- -> Cu','2H2O + 2e- -> H2 + 2OH-','2H2O -> 2H2 + O2','CuO + H2SO4 -> CuSO4 + H2O','Cu(OH)2 + H2SO4 -> CuSO4 + 2H2O']:
    assert check_equation(e) is None, e
ne = 1.00*3860/9.65e4; assert abs(ne-0.0400)<1e-12
nCu = 0.100*0.200; assert abs(nCu-0.0200)<1e-12
m1 = nCu*(63.5+16); chk('m1',m1,1.59,3)
m2 = m1 + (ne/2)*18.0; chk('m2',m2,1.95,3)
# atoms lost after stage 2: Cu 0.02, O 0.02+0.02, H 0.04
lostCu, lostO, lostH = nCu, nCu + ne/2, 2*(ne/2)
assert abs(lostO-0.04)<1e-12 and abs(lostH-0.04)<1e-12
chk('Cu(OH)2',0.0200*(63.5+2*17.0),1.95,3)

print('ALL OK')
