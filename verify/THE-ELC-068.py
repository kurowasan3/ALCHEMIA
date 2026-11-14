ITEM = '/home/claude/chemx/tokodai/out/items/THE-ELC-068/item.json'
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
x,h=sp.symbols('x h')
s=sp.solve([x+0.100+h/2-0.480, 24*x+5.90+1.0*h-11.06],[x,h]); print(s)
chk('Mg',float(s[x]),0.200,3); chk('H',float(s[h]),0.360,3)
chk('H%',float(s[h])/11.06*100,3.3,2)
chk('Q',float(s[h])*96500,34700,3); chk('ratio',0.360/0.960,0.38,2)
chk('Ni',5.90/59,0.100,3)
print('ALL OK')
