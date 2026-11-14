ITEM = '/home/claude/chem6/out/items/INO-NIT-019/item.json'
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
x=sp.symbols('x')
nMg=sp.Rational(120,100)/24
sol=sp.solve(sp.Eq(40*(nMg-x)+sp.Rational(100,3)*x, sp.Rational(190,100)),x)[0]
print('nMg',nMg,'x',sol,'%',sol/nMg*100)
assert sol==sp.Rational(15,1000) and sol/nMg*100==30
assert sig(float(nMg*40),3)==2.00 and sig(float(nMg/3*100),3)==1.67
nNH3=sol/3*2; Hres=sp.Rational(2,100)-nNH3; vol=Hres/sp.Rational(2,10)*1000
print('NH3',nNH3,'NaOH mL',vol); assert vol==50

print('全検証OK:', ITEM)
