ITEM = '/home/claude/chem4/outh/items/THE-EQL-070/item.json'
# 検証スクリプト共通部（グループ g04）  python3 -I この.py
import sys, json, re
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
# 小問数と答えの番号の一致
nprob = it['problem'].count('\\item')
nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
if nprob:
    assert [int(n) for n in nums] == list(range(1, nprob + 1)), (nprob, nums)
print('小問数OK:', nprob)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import math
k=0.693/50/0.10**2; print('k',k); assert sig(k,2)==1.4
assert round(math.log2((50/12.5)),6)==2
def half(A0,B0,k=k,h=0.002):
    A,B,t=A0,B0,0.0
    f=lambda A,B: -k*A*B*B
    while A>A0/2:
        k1=f(A,B); k2=f(A+h/2*k1,B+h*k1); k3=f(A+h/2*k2,B+h*k2); k4=f(A+h*k3,B+2*h*k3)
        dA=h/6*(k1+2*k2+2*k3+k4); A+=dA; B+=2*dA; t+=h
    return t
hs=[half(1e-3,0.10),half(1e-3,0.20),half(5e-4,0.10)]; print(hs)
assert abs(hs[0]-50)<0.6 and abs(hs[1]-12.5)<0.1 and abs(hs[2]-50)<0.4
r=2/(1-math.sqrt(0.95)); print('ratio',r, r*1e-3); assert sig(r*1e-3,2)==0.079
assert all(x>=r for x in (100,200,200))
v0=k*0.010*0.020**2; print('v0',v0); assert sig(v0,2)==5.5e-6
# 化学量論比：半減期
th=half(0.010,0.020,h=0.05); th2=half(0.020,0.040,h=0.02); print(th, 3/(8*k*0.01**2), th2/th)
assert abs(th-3/(8*k*0.01**2))/th<0.002 and abs(th2/th-0.25)<0.002
print('全検証OK:', ITEM)
