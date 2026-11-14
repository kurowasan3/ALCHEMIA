ITEM = '/home/claude/chem6/outh/items/INO-NIT-022/item.json'
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
n=sp.symbols('n')
Na=4+(n-2); H=2+2*(n-2)-2*(n-1); Ox=8+4*(n-2)-(n-1); P=n
assert sp.simplify(Na-(n+2))==0 and sp.simplify(H)==0 and sp.simplify(Ox-(3*n+1))==0
for e in ['2Na2HPO4 + NaH2PO4 -> Na5P3O10 + 2H2O','2Na2HPO4 -> Na4P2O7 + H2O','6Na2HPO4 + 3NaH2PO4 -> 3Na5P3O10 + 6H2O']:
    assert check_equation(e) is None
x,y=sp.symbols('x y')
s=sp.solve([142*x+120*y-sp.Rational(644,100), x+2*y-sp.Rational(8,100)],[x,y]); print(s)
assert s[x]==sp.Rational(2,100) and s[y]==sp.Rational(3,100)
assert 18*(s[x]+2*s[y])/2==sp.Rational(72,100)
Pn=s[x]+s[y]; chains=s[x]/2; print('avg',Pn/chains); assert Pn/chains==5
assert (2*s[x]+s[y])/Pn==sp.Rational(7,5)
L=sp.Rational(630,1000)/18; print('L',L)
OHres=(s[x]+2*s[y])-2*L; ch2=Pn-L; avg2=Pn/ch2; print('OH',OHres,'avg',avg2,float(avg2))
assert OHres==sp.Rational(1,100) and avg2==sp.Rational(10,3) and sig(float(avg2),2)==3.3
r=sp.solve(sp.Eq((2*x+y)/(x+y),sp.Rational(11,10)),y)[0]; print('y/x',r/x); assert r/x==9
assert sp.Rational(20+2,20)==sp.Rational(11,10)

print('全検証OK:', ITEM)
