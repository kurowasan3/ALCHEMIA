ITEM = '/home/claude/chem4/out/items/THE-ION-065/item.json'
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
K=sp.Rational(54,10)**2; print('K',float(K)); assert sig(float(K),2)==29
s=sp.symbols('s')
sol=[r for r in sp.solve((sp.Rational(54,10)-s)*(sp.Rational(64,10)-s)-K,s) if 0<r<sp.Rational(54,10)]
print(sol,[float(r) for r in sol]); assert len(sol)==1
m=float(sol[0])*58.5; print('mass',m); assert sig(m,2)==28
# (3) 0.10 mol ずつ加えたときのイオン積
c=5.4
cases={'a':(c,c),'b':(c+0.1,c),'c':(c,c+0.1),'d':(c/1.1,c/1.1)}
ppt=[k for k,(na,cl) in cases.items() if na*cl>float(K)+1e-12]; print(ppt); assert ppt==['b','c']
print('全検証OK:', ITEM)
