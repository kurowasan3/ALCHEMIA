import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-022/item.json'
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

for e in ['LiC6 -> C6 + Li+ + e-', 'CoO2 + Li+ + e- -> LiCoO2', '2Li + 2H2O -> 2LiOH + H2']:
    assert check_equation(e) is None
close(Q(1, 10)*Q(69, 10), 0.69)
r = (Q(69, 10) + 72)/Q(69, 10); print('ratio', float(r)); assert sig(r, 2) == 11
# oxidation numbers: CoO2 -> +4 ; LiCoO2 -> +3
assert 0 - 2*(-2) == 4 and 0 - 1 - 2*(-2) == 3
print('THE-ELC-022 OK')
