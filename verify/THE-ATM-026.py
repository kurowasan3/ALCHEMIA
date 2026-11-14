ITEM = '/home/claude/chemx/x1/out/items/THE-ATM-026/item.json'
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

rC=1.1/98.9; rH=0.015/99.985; r17=0.04/99.76; r18=0.20/99.76
a = {'13C':2*rC,'2H':6*rH,'17O':r17}
assert max(a,key=a.get)=='13C'; chk('M+1 13C %',2*rC*100,2.2,2)
b = {'18O':r18,'13C2':rC**2,'13C2H':12*rC*rH,'17O2H':6*r17*rH,'13C17O':2*rC*r17,'2H2':15*rH**2}
assert max(b,key=b.get)=='18O'; chk('18O %',r18*100,0.20,2)
chk('13C2 %',rC**2*100,0.012,2)
chk('n',0.067/rC,6.0,2)

print('ALL OK')
