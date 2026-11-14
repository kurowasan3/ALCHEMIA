ITEM = '/home/claude/chem5/outh/items/THE-ION-069/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from sympy import Rational as R, nsimplify
Ka = R(3,10**4); Hs = R(1,100); Hp = R(4,10**8)
r1 = (1+Ka/Hp)/(1+Ka/Hs); print(float(r1)); assert sig(float(r1),2)==7300
Kb = R(1,10**9)
r2 = (1+Hs/Kb)/(1+Hp/Kb); print(float(r2)); assert sig(float(r2),2)==2.4e5
u1 = (1+Ka/R(1,10**5))/(1+Ka/Hp); u2 = (1+Ka/R(1,10**8))/(1+Ka/Hp)
print(float(u1), float(u2), float(u2/u1)); assert sig(float(u1),2)==4.1e-3 and sig(float(u2),2)==4.0 and 960<float(u2/u1)<980
approx = (1+Ka/Hp)/1; dev = float(approx/r1-1); print('ずれ', dev); assert abs(dev-0.03)<1e-9 and sig(float(approx),2)==7500
Ka2 = 10*Ka; r4 = (1+Ka2/Hp)/(1+Ka2/Hs); print(float(r4), float(r4/r1)); assert sig(float(r4),3)==5.77e4 and sig(float(r4/r1),2)==7.9
print('血漿側の1を無視した誤差', float(1/(Ka/Hp)))
print('ALL OK')
