ITEM = '/home/claude/chem2/outh/items/THE-ION-051/item.json'
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

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
check_all_equations()
assert check_equation('AgCl(s) <=> Ag+ + Cl-') is None and check_equation('Ag+ + 2Cl- <=> [AgCl2]-') is None
from sympy import Rational as R, symbols, sqrt, diff, solve
x = symbols('x', positive=True)
K = R(18, 10**11); b = 10**5
S = K/x + K*b*x
s0 = sqrt(K); print('(2)', float(s0), float(K*b*s0)); assert sig(float(s0),2)==1.3e-5 and sig(float(K*b*s0),2)==2.4e-10
xm = solve(diff(S, x), x)[0]; Sm = S.subs(x, xm)
print('(3)', float(xm), float(Sm)); assert sig(float(xm),2)==3.2e-3 and sig(float(Sm),2)==1.1e-7
S1 = S.subs(x, 1); print('(4)', float(S1), float(K*b/S1)); assert sig(float(S1),2)==1.8e-5 and float(S1) > float(s0)
r = sorted(float(v) for v in solve(S - R(2, 10**6), x)); print('(5) exact', r)
assert sig(r[0],2)==9.0e-5 and sig(r[1],2)==0.11
assert sig(float(K/R(2,10**6)),2)==9.0e-5 and sig(float(R(2,10**6)/(K*b)),2)==0.11

print('全検証OK:', ITEM)
