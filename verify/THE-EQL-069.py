ITEM = '/home/claude/chem4/outh/items/THE-EQL-069/item.json'
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
kd=math.log(2)/5; k0=math.log(4)*kd
print('k0',k0); assert sig(kd,2)==0.14 and sig(k0,2)==0.19
nd=math.exp(-k0*5); print('no deact',nd); assert sig(nd,2)==0.38 and abs(math.exp(-0.96)-0.38)<0.005
def run(adds, C0=1.0, T=200, h=0.001, probe=()):
    # adds: list of (time, k0_i)
    C=C0; t=0.0; out={}
    def k(t): return sum(a*math.exp(-kd*(t-ti)) for ti,a in adds if t>=ti-1e-12)
    n=int(T/h)
    for i in range(n):
        for p in probe:
            if abs(t-p)<h/2: out[p]=C
        f=lambda tt,cc: -k(tt)*cc
        k1=f(t,C); k2=f(t+h/2,C+h/2*k1); k3=f(t+h/2,C+h/2*k2); k4=f(t+h,C+h*k3)
        C+=h/6*(k1+2*k2+2*k3+k4); t+=h
    return C,out
Cinf,o=run([(0,k0)],probe=(5.0,)); print(Cinf,o)
assert abs(Cinf-0.25)<1e-4 and abs(o[5.0]-0.5)<2e-3
# 式(*)
t=7.0; assert abs(math.exp(-(k0/kd)*(1-math.exp(-kd*t))) - run([(0,k0)],T=t)[0])<1e-4
m=math.log(100)/math.log(4); print('m',m); assert sig(m,2)==3.3
Cinf2,o2=run([(0,k0/2),(10,k0/2)],probe=(9.999,)); print(Cinf2,o2)
assert abs(Cinf2-0.25)<1e-4
c10=2**-0.75; assert sig(c10,2)==0.59 or sig(1/1.68,2)==0.60
assert abs(list(o2.values())[0]-c10)<2e-3
Cd,_=run([(0,k0)],C0=0.5); assert abs(Cd/0.5-0.25)<1e-4
print('全検証OK:', ITEM)
