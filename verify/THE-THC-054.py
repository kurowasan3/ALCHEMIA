ITEM = '/home/claude/chemx/tokodai/out/items/THE-THC-054/item.json'
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
cP,cE,cY,cH=2220,2058,1938,286
fin=3.00*7/15; used=3.00-fin; chk('used',used,1.60,3)
b,c=sp.symbols('b c'); s=sp.solve([b+c-1,b+2*c-used],[b,c]); print(s)
QX=cP+(2.00-used)*cH; QY=float(s[b])*cE+float(s[c])*cY
chk('QX',QX,2330,3); chk('QY',QY,1990,3)
hyd=QY+2*cH-QX; chk('hyd',hyd,224,3)
assert abs(hyd-(float(s[b])*(cE+cH-cP)+float(s[c])*(cY+2*cH-cP)))<1e-9
co2={'P':cP/3,'E':cE/3,'Y':cY/3}; o2={'P':cP/5,'E':cE/4.5,'Y':cY/4}
assert max(co2,key=co2.get)=='P' and max(o2,key=o2.get)=='Y'; print(co2,o2)
print('ALL OK')
