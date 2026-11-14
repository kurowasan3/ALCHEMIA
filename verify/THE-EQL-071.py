ITEM = '/home/claude/chem4/outh/items/THE-EQL-071/item.json'
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
import sympy as sp, math
x,P,K=sp.symbols('x P K',positive=True)
Kp=(x*P)**2/(((1-x)/4*P)*(3*(1-x)/4*P)**3)
assert sp.simplify(Kp-256*x**2/(27*(1-x)**4*P**2))==0
lhs=sp.sqrt(Kp)*3*sp.sqrt(3)/16*P
assert sp.simplify(lhs - x/(1-x)**2)==0
c=lambda v: v/(1-v)**2
A=3*math.sqrt(3)/16
K773=(c(0.100)/(A*1e7))**2; print('K773',K773); assert sig(K773,2)==1.4e-15
def solve(cc):
    r=[float(s) for s in sp.solve(sp.Eq(x/(1-x)**2,cc),x) if 0<float(s)<1]; assert len(r)==1; return r[0]
x3=solve(c(0.1)*3); print('x3',x3); assert sig(x3,2)==0.22
lnr=92e3/8.31*(1/673-1/773); print('ln ratio',lnr); assert abs(lnr-2.13)<0.005
x400=solve(c(0.1)*math.sqrt(math.exp(lnr))); print('x400',x400); assert sig(x400,2)==0.22
x6=solve(c(0.1)/10); print('x 1e6',x6); assert sig(x6,2)==0.012
print('sqrtK*P ratio', math.sqrt(math.exp(lnr)), 3)
print('全検証OK:', ITEM)
