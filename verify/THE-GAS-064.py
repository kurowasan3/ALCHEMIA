ITEM = '/home/claude/chemx/x1/out/items/THE-GAS-064/item.json'
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

R=8310;T=300;P0=1.00e5;pw=3.6e3
nHe=(P0-pw)*4.80/(R*T); nX=(P0-pw)*1.20/(R*T)
chk('nHe',nHe,0.186,3); chk('nX',nX,0.0464,3)
VHe=nHe*R*T/1.03e5; VX=nX*R*T/1.02e5
chk('VHe',VHe,4.49,3); chk('VX',VX,1.13,3)
rho=P0*29.0/(R*T); air=(VHe+VX)*rho; chk('air',air,6.55,3)
mX=air-2*1.75-nHe*4.0; chk('mX',mX,2.30,3)
M=mX/nX; chk('M',M,50,2)
x=(58-M)/14; print('x propane',x); assert 0.5<x<1

print('ALL OK')
