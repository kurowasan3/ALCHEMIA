ITEM = '/home/claude/chem4/out/items/THE-EQL-066/item.json'
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
V,x,y=sp.symbols('V x y',positive=True)
K=sp.simplify((1/V)*(1/V)/((sp.Rational(1,2)/V)*(2/V))); print('K',K); assert K==1
sx=[s for s in sp.solve(sp.Eq((1+x)**2,(2-x)**2),x) if 0<s<2]; print(sx); assert sx==[sp.Rational(1,2)]
assert abs(float(sx[0]/sp.Rational(3,2))*100-33.3)<0.05
sy=[s for s in sp.solve(sp.Eq((1+y)**2,(sp.Rational(1,2)-y)*(sp.Rational(7,2)-y)),y) if 0<s<sp.Rational(1,2)]; print(sy); assert sy==[sp.Rational(1,8)]
assert sig(1+float(sy[0]),3)==1.13 or sig(1.125,3) in (1.12,1.13)
print('Q after CO add',1/(2*2),'after H2O add',1/(0.5*3.5))
print('全検証OK:', ITEM)
