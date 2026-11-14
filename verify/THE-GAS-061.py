ITEM = '/home/claude/chemx/tokodai/outh/items/THE-GAS-061/item.json'
# 検証スクリプト（東工大 g03）  python3 -I この.py
# 共通の検証ヘッダ（各検証スクリプトにそのまま埋め込む）
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
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')
import scipy.optimize as so
ps={275:7.0,284:13.0,285:13.9,286:14.8,287:15.8,288:16.9,289:18.0,290:19.2}
f=lambda K:6.1078*10**(7.5*(K-273.15)/(K-273.15+237.3))
for T,p in ps.items(): assert abs(round(f(T),1)-p)<0.051,(T,f(T))
chk('Tdry',300-25.0,275,3); chk('pw',750*0.03,22.5,3); assert 22.5>ps[275]
n=lambda p,nd=97.0,PC=750: nd*p/(PC-p)
R=lambda T,w0=3.0,nd=97.0,PC=750: 300-(1000-PC)/10+44/3*(w0-n(ps[T],nd,PC))
for T in (286,287,288,289): print(T, round(n(ps[T]),3), round(R(T),1), round(R(T)-T,2))
assert R(287)-287>0>R(288)-288
g1,g2=R(287)-287,R(288)-288; TC=287+g1/(g1-g2); chk('TC',TC,287,3); chk('TD',TC+25,312,3)
def solve(w0,nd,PC):
    g=lambda T:300-(1000-PC)/10+44/3*(w0-nd*f(T)/(PC-f(T)))-T
    T=so.brentq(g,250,300); return T,w0-nd*f(T)/(PC-f(T)),T+(1000-PC)/10
b=solve(3.0,97.0,750); a=solve(3.0,97.0,700); i=solve(2.5,97.5,750); print(b,a,i)
assert abs(b[0]-TC)<0.05 and a[2]>b[2] and i[2]<b[2]
chk('dT',44/3*b[1],12,2)
print('ALL OK')
