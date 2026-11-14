ITEM = '/home/claude/chemx/x1/outh/items/THE-EQL-142/item.json'
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
ct=0.0200
def solve_nh3(c):  # (4-y)/(10+y)*ct=c
    return (4-10*c/ct)/(1+c/ct)
def solve_h2(c):   # 1.5y/(10+y)*ct=c
    return 10*c/ct/(1.5-c/ct)
y1=solve_nh3(5.00e-3); chk('y1%',y1/4*100,30,2)
t1=10+y1
chk('N2',y1/2/t1*ct,1.07e-3,3); chk('H2',1.5*y1/t1*ct,3.21e-3,3); chk('Ar',6/t1*ct,1.07e-2,3)
y2=solve_h2(5.00e-3); chk('y2%',y2/4*100,50,2)
t2=10+y2
Kc=(y2/2/t2*ct)*(1.5*y2/t2*ct)**3/((4-y2)/t2*ct)**2; chk('Kc',Kc,1.88e-5,3)
Q1=(y1/2/t1*ct)*(1.5*y1/t1*ct)**3/((4-y1)/t1*ct)**2; chk('Q1',Q1,1.42e-6,3); assert Q1<Kc
chk('体積比',t2/10,1.2,2); chk('エ',(1-5/8)*100,37.5,3)
chk('Ar出口',6/t2*ct,1.00e-2,3)

print('ALL OK')
