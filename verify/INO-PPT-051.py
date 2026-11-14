ITEM = '/home/claude/chemx/x1/outh/items/INO-PPT-051/item.json'
# 検証スクリプト（x1 g04）  python3 -I この.py
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
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

import itertools
perms=[]
for r in range(4):
    perms.append(lambda a,r=r: a[r:]+a[:r])
    perms.append(lambda a,r=r: (a[r:]+a[:r])[::-1])
def canon(a): return min(p(a) for p in perms)
def inv(a): return tuple((e, 'D' if c=='L' else 'L') for e,c in a)
def build(chirs):
    out=set()
    for c1,c2 in chirs:
        for o1,o2 in itertools.product(range(2),repeat=2):
            l1=[('N',c1),('O',c1)]; l2=[('N',c2),('O',c2)]
            if o1: l1=l1[::-1]
            if o2: l2=l2[::-1]
            out.add(canon(tuple(l1+l2)))
    return out
onlyL=build([('L','L')]); assert len(onlyL)==2
allc=build([(x,y) for x in 'LD' for y in 'LD']); assert len(allc)==6, len(allc)
achiral=[a for a in allc if canon(inv(a))==a]
assert len(achiral)==2 and all({c for _,c in a}=={'L','D'} for a in achiral)
# L だけの2種は互いに鏡像でない
a,b=sorted(onlyL); assert canon(inv(a))!=b and canon(inv(a)) not in onlyL
print('立体異性体 6, 光学不活性 2 OK')
for p,exp in [(0.5,50),(0.75,37.5)]:
    chk('不活性%', 2*p*(1-p)*100, exp, 3)

print('ALL OK')
