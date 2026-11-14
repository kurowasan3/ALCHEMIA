ITEM = '/home/claude/chemx/x1/out/items/THE-EQL-136/item.json'
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

assert check_equation('2SO3 -> 2SO2 + O2') is None
Kp=2.0e-3; P=1.0e5
chk("Kp'",Kp**2,4.0e-6,2)
def bis(f,lo,hi):
    for _ in range(300):
        m=(lo+hi)/2
        if f(lo)*f(m)<=0: hi=m
        else: lo=m
    return m
# (2) 純 SO3：解離量 2x（SO2），x（O2），全圧一定で分圧はモル分率×P
def g2(x):  # 物質量 SO3 1-2x, SO2 2x, O2 x
    t=1+x; return (2*x/t*P)*((x/t*P)**0.5)/((1-2*x)/t*P)-Kp
x=bis(g2,1e-15,0.4); pO2=x/(1+x)*P; chk('pO2(2)',pO2,22,2)
chk('解離%',2*x*100,0.043,2)
def g3(x):  # SO3 0.5-2x, SO2 0.5+2x, O2 x
    t=1+x; return ((0.5+2*x)/t*P)*((x/t*P)**0.5)/((0.5-2*x)/t*P)-Kp
y=bis(g3,1e-30,0.1); chk('pO2(3)',y/(1+y)*P,4.0e-6,2)

print('ALL OK')
