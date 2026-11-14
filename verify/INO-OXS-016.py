ITEM = '/home/claude/chem6/outh/items/INO-OXS-016/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
print('検査した式の数:', n_eq)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


import sympy as sp
a=sp.symbols('a')
def eqsolve(so2,o2,n2,KP):
    tot=so2+o2+n2
    f=(a/(1-a))**2-KP*(o2-so2*a/2)/(tot-so2*a/2)
    return float(sp.nsolve(f,a,(0.5,0.99999),solver='bisect'))
KP=0.16*1.0e5
a1=eqsolve(10,11,79,KP); so2=10*(1-a1); o2=11-5*a1
a2=eqsolve(so2,o2,79,KP); tot=1-(1-a1)*(1-a2); aP=eqsolve(10,11,79,2*KP)
print(a1,a2,tot,aP)
assert sig(a1,3)==0.970 and sig(a2,3)==0.971 and sig(tot,3)==0.999 and sig(aP,3)==0.978
# 近似のくり返し
import math
r=math.sqrt(KP*6/95); al=r/(1+r)
for _ in range(3):
    xo=(11-5*al)/(100-5*al); r=math.sqrt(KP*xo); al=r/(1+r)
print('iter',al); assert sig(al,3)==0.970
ratio=(1-aP)/((1-a1)*(1-a2)); print('ratio',ratio); assert 24<ratio<26
assert check_equation('2SO2 + O2 <=> 2SO3') is None

print('全検証OK:', ITEM)
