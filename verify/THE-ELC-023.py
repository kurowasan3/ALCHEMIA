import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-023/item.json'
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

for e in ['VO2+ + 2H+ + e- -> VO^2+ + H2O', 'V^2+ -> V^3+ + e-', 'VO2+ + V^2+ + 2H+ -> VO^2+ + V^3+ + H2O',
          'VO^2+ + H2O -> VO2+ + 2H+ + e-', 'V^3+ + e- -> V^2+']:
    assert check_equation(e) is None
x = symbols('x'); assert solve(x + 2*(-2) - 1, x) == [5]; assert solve(x - 2 - 2, x) == [4]
n = (Q(14, 10) - Q(2, 10))*2; close(n, 2.4)
q = n*96500; print(float(q)); assert sig(q, 2) == 2.3e5
# charge change per e- : positive side cations VO2+(+1)->VO^2+(+2), -2H+  => -1 ; negative side +1
assert (2 - 1 - 2) == -1
print('THE-ELC-023 OK')
