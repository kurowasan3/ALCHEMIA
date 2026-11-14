ITEM = '/home/claude/chem4/outh/items/THE-EQL-072/item.json'
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
v,K=sp.symbols('v K',positive=True)
def r_of(Kc,u):
    S=u+2*Kc*u*u
    sol=[s for s in sp.solve(sp.Eq(v+2*Kc*v**2,2*S),v) if s>0]; return float(sol[0])/float(u)
r1=r_of(sp.Integer(200),sp.Rational(1,100)); print('r1',r1); assert sig(r1,3)==1.46
prev=None
for e in range(-4,5):
    Ku=10.0**e; r=r_of(sp.Float(Ku/0.01),sp.Rational(1,100))
    assert math.sqrt(2)<r<2; 
    if prev: assert r<prev
    prev=r; print(Ku, r)
assert r_of(sp.Float(1e6),sp.Rational(1,100))-math.sqrt(2)<1e-3
Ks=sp.solve(sp.Eq(sp.Rational(15,1000)+2*K*sp.Rational(15,1000)**2, 2*(sp.Rational(1,100)+2*K*sp.Rational(1,100)**2)),K); print(Ks); assert Ks==[100]
assert r1/2<1
print('全検証OK:', ITEM)
