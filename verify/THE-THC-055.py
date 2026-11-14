ITEM = '/home/claude/chemx/tokodai/out/items/THE-THC-055/item.json'
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
import sympy as sp
u=1.5; chk('CO',2*u,3.00,3); chk('CO2',u,1.50,3)
Qa=2*u*(394-283)+u*394; chk('Qa',Qa,924,3)
y,r=sp.symbols('y r')
ys=sp.solve(3.00+y-3.00*(3*u),y)[0]; chk('y',float(ys),10.5,3)
rs=sp.solve(3.00*283+394*r-2425,r)[0]; chk('r',float(rs),4.00,3)
tot=3*u+float(rs); chk('tot',tot,8.50,3)
assert float(ys)-1.5-float(rs)>0
chk('sum',Qa+2425,3350,3); assert abs(Qa+2425-tot*394)<1e-6
print('ALL OK')
