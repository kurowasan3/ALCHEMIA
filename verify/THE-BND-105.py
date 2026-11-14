# 検証スクリプト THE-BND-105   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-BND-105/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from math import sqrt
from itertools import product
def svals(rule, n):
    S = set()
    for h, k, l in product(range(0, 6), repeat=3):
        if (h, k, l) == (0, 0, 0): continue
        if rule(h, k, l): S.add(h*h+k*k+l*l)
    return sorted(S)[:n]
sc = svals(lambda h,k,l: True, 7); bcc = svals(lambda h,k,l: (h+k+l) % 2 == 0, 7); fcc = svals(lambda h,k,l: len({h%2,k%2,l%2}) == 1, 6)
print(sc, bcc, fcc); assert sc == [1,2,3,4,5,6,8] and bcc[:6] == [2,4,6,8,10,12] and bcc[6] == 14 and fcc == [3,4,8,11,12,16]
lam = 154; a_true = 314.7
obs = [round(lam**2*s/(4*a_true**2), 4) for s in [2,4,6,8,10,12]]
print(obs); assert obs == [0.1197, 0.2395, 0.3592, 0.4789, 0.5987, 0.7184]
ratio = [o/obs[0] for o in obs]; assert all(abs(r - (i+1)) < 0.01 for i, r in enumerate(ratio))
assert abs(fcc[1]/fcc[0] - 2.0) > 0.5
s1 = sin1 = sqrt(obs[0])
a_sc = lam/(2*sin1); a_bcc = lam*sqrt(2)/(2*sin1)
Z = lambda a: 10.2*(a*1e-10)**3*6.02e23/95.9
print(a_sc, Z(a_sc), a_bcc, Z(a_bcc)); assert sig(a_sc,3) == 223 and sig(Z(a_sc),2) == 0.71 and sig(a_bcc,3) == 315 and sig(Z(a_bcc),2) == 2.0
r = sqrt(3)*a_bcc/4; print(r); assert sig(r,3) == 136
s7b = obs[0]*14/2; s7s = obs[0]*8; print(s7b, s7s); assert sig(s7b,3) == 0.838 and sig(s7s,3) == 0.958
# 7 は3つの平方数の和で表せない
assert 7 not in {h*h+k*k+l*l for h, k, l in product(range(3), repeat=3)}

print('検証OK')
