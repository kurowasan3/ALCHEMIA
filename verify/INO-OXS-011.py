import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/INO-OXS-011/item.json'
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

h2so4 = Q(500)*Q(98, 100)/98*1000/1000  # kmol? use mol
n = Q(50000)*Q(98, 100)/98; close(n, 500)
close(n/2*120/1000, 30.0)
r = (n/2*120)/(n*32); print('ratio', float(r)); assert sig(r, 2) == 1.9
a = symbols('a'); av = solve(100 - a - 96, a)[0]
conv = 2*av/10*100; close(conv, 80); assert av <= 10
for e in ['2SO2 + O2 -> 2SO3', 'SO3 + H2O -> H2SO4']:
    assert check_equation(e) is None
print('INO-OXS-011 OK')
