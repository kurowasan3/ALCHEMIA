ITEM = '/home/claude/chemx/x1/out/items/THE-EQL-135/item.json'
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

assert check_equation('2NH3 -> N2 + 3H2') is None
x=2.70*2/3; NH3=3.00-x; N2=0.900; tot=NH3+N2+2.70+5.00
chk('NH3',NH3,1.20,3); chk('tot',tot,9.80,3)
Q=400*tot/8.00; chk('Q',Q,490,3)
Kc=(N2/Q)*(2.70/Q)**3/(NH3/Q)**2; chk('Kc',Kc,5.12e-5,3)
ctot=8.00/400
def f(xx,Ar,F0=3.00):
    t=F0+xx+Ar
    return (xx/2/t)*(1.5*xx/t)**3/((F0-xx)/t)**2*ctot**2-Kc
def solve(Ar,F0=3.00):
    lo,hi=1e-9,F0-1e-12
    for _ in range(200):
        m=(lo+hi)/2
        if f(m,Ar,F0)>0: hi=m
        else: lo=m
    return m
chk('H2(Ar=5)',1.5*solve(5.00),2.70,3)
chk('H2(Ar=15)',1.5*solve(15.0),3.24,3)
chk('H2(2倍)',1.5*solve(10.0,6.00),5.40,3)

print('ALL OK')
