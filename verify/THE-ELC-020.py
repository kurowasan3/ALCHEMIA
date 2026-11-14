import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-020/item.json'
d = json.load(open(ITEM, encoding='utf-8'))


def check_all_ce(d):
    txt = d['problem'] + '\n' + d['answer'] + '\n' + d['solution']
    n = 0
    for m in re.finditer(r'\\ce\{([^{}]*)\}', txt):
        s = m.group(1)
        arrows = re.findall(r'<=>|->', s)
        if len(arrows) != 1 or 'boxed' in s or re.search(r'\bR\b|\bn\b', s):
            continue
        s = re.sub(r'M(?![a-z])', 'Xe', s)  # 一般の水素吸蔵合金 M を仮の元素で置き換えて検査
        r = check_equation(s)
        assert r is None, (s, r)
        n += 1
    print('equations checked:', n)


def close(a, b, rel=5e-3):
    a = float(a); b = float(b)
    assert abs(a - b) <= rel * max(abs(b), 1e-12), (a, b)


def sig(x, k):
    from math import log10, floor
    x = float(x)
    return round(x, -int(floor(log10(abs(x)))) + (k - 1))


check_all_ce(d)

gas = Q(168, 1000)/Q(224, 10); ne = gas*Q(4, 3); close(ne, 0.01)
t = ne*96500/Q(4, 10)/60; print('min', float(t)); assert sig(t, 3) == 40.2
close(ne*98, 0.980)
cu_after = Q(5, 100)*Q(12, 1000)/Q(10, 1000); close(cu_after, 0.06)
dec = ne/2/Q(2, 10); c0 = cu_after + dec; close(c0, 0.085)
assert cu_after > 0
for e in ['Pb + PbO2 + 2H2SO4 -> 2PbSO4 + 2H2O', '2H+ + 2e- -> H2']:
    assert check_equation(e) is None
print('THE-ELC-020 OK')
