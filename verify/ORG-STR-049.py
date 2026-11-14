ITEM = '/home/claude/chem6/out/items/ORG-STR-049/item.json'
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


import itertools, sympy as sp
G={'CH4':(1,4),'C2H4':(2,4),'C2H6':(2,6),'C3H6':(3,6),'C3H8':(3,8),'C4H10':(4,10)}
p1=[(a,b) for a,b in itertools.combinations(G,2) if G[a][0]==G[b][0]==3]
print(p1); assert p1==[('C3H6','C3H8')]
sameC=[(a,b) for a,b in itertools.combinations(G,2) if G[a][0]==G[b][0]]; print('sameC',sameC)
q=[(a,b) for a,b in itertools.combinations(G,2) if G[a][1]==G[b][1] and G[a][0]!=G[b][0]]
print(q); assert q==[('CH4','C2H4'),('C2H6','C3H6')]
a,b=sp.symbols('a b'); s=sp.solve([a+b-10,3*a+4*b-36],[a,b]); print(s); assert s=={a:4,b:6}
o2=s[a]*sp.Rational(9,2)+s[b]*5; assert o2==48 and (30*2+36)==2*o2
for e in ['2C3H6 + 9O2 -> 6CO2 + 6H2O','C3H8 + 5O2 -> 3CO2 + 4H2O','C3H6 + H2 -> C3H8']:
    assert check_equation(e) is None

print('全検証OK:', ITEM)
