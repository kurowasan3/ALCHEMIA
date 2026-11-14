ITEM = '/home/claude/chemx/x1/out/items/THE-BND-090/item.json'
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

a=0.409e-7; M=108; NA=6.02e23; s3=1.73
rho=4*M/(NA*a**3); chk('密度', rho, 10.5, 3)
t=2.16/rho/2.00e4; chk('厚さnm', t*1e7, 103, 3)
d=a/s3; chk('間隔', d*1e7, 0.236, 3)
chk('層数', t/d, 436, 3)
sg=4/(s3*a*a); chk("面密度", sg, 1.38e15, 3)
N=2.16/M*NA; chk('層数2', N/(2.00e4*sg), 436, 3)
assert abs(t/d - N/(2.00e4*sg)) < 1e-6
# 六方最密: 面密度・層間隔は同じ（2r=a/√2 の最密層，間隔 sqrt(2/3)*2r）
r2=a/2**0.5; assert abs((2/3)**0.5*r2 - a/3**0.5) < 1e-20

print('ALL OK')
