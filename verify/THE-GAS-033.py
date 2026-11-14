ITEM = '/home/claude/chem2/out/items/THE-GAS-033/item.json'
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
assert not check_all_equations()
from sympy import Rational as R, symbols, solve, diff
rc, ra, rar = R('1.2505'), R('1.2569'), R('1.7837')
x = (ra - rc) / (rar - rc); print('(2)', float(x)); assert sig(float(x)*100, 3) == 1.20
air = R('0.790') * x; print('(3)', float(air)); assert sig(float(air)*100, 3) == 0.948
rho = symbols('rho'); xx = (ra - rho) / (rar - rho)
assert diff(xx, rho).subs(rho, rc) < 0
x_bias = (ra - (rc*R('0.999') + R('1.4290')*R('0.001'))) / (rar - rc); assert x_bias < x
assert R('0.1785') < rc
for g, M in [(rc, 28.01), (rar, 39.95), (R('1.4290'), 32.0), (R('0.1785'), 4.00)]:
    assert abs(float(g)*22.41 - M) / M < 0.003, (g, M)

print('全検証OK:', ITEM)
