ITEM = '/home/claude/chemx/tokodai/out/items/THE-GAS-059/item.json'
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
f=lambda K:6.1078*10**(7.5*(K-273.15)/(K-273.15+237.3))
tab={0:(1000,300,35.3),200:(980,298,31.4),400:(960,296,27.8),600:(940,294,24.6)}
for z,(P,T,ps) in tab.items(): assert abs(round(f(T),1)-ps)<0.051, (T,f(T))
assert abs(round(f(284),1)-13.0)<0.051
x=0.0300; chk('p0',1000*x,30.0,3); chk('p400',960*x,28.8,3)
d200=980*x-31.4; d400=960*x-27.8; assert d200<0<d400
chk('z',200+200*(-d200)/(d400-d200),330,2)
n=97.0*13.0/(780-13.0); chk('n',n,1.64,3); chk('cond',3.00-n,1.36,3)
chk('pD',1000*n/(97.0+n),16.7,3)
assert all(P*0.02<ps for z,(P,T,ps) in tab.items())
print('ALL OK')
