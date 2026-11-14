ITEM = '/home/claude/chem4/out/items/THE-EQL-068/item.json'
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
import sympy as sp
x,y,z,V=sp.symbols('x y z V',positive=True)
xs=sp.solve(sp.Eq((12-2*x)/12,sp.Rational(3,4)),x)[0]; assert xs==sp.Rational(3,2)
n2,h2,nh3=4-xs,8-3*xs,2*xs; print(n2,h2,nh3)
assert sig(float(n2/h2),2)==0.71
Kc=sp.simplify((nh3/V)**2/((n2/V)*(h2/V)**3)); print(Kc); assert sp.simplify(Kc-sp.Rational(144,1715)*V**2)==0
zs=sp.solve(sp.Eq(h2-3*z,2*z),z)[0]; assert zs==sp.Rational(7,10)
ys=sp.solve(sp.Eq((2*zs/V)**2/(((n2+y-zs)/V)*((h2-3*zs)/V)**3),Kc),y); print(ys,[float(v) for v in ys])
assert len(ys)==1 and sig(float(ys[0]),2)==6.7
print('全検証OK:', ITEM)
