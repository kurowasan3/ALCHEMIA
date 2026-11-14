import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-EQL-032/item.json'
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

K = Q(6, 100)/(Q(4, 100)*Q(1, 2)); close(K, 3.0)
cases = {'2': (Q(8, 100), Q(2, 100), Q(1, 2)), '3': (Q(2, 100), Q(8, 100), Q(1, 2)),
         '4': (Q(4, 100), Q(6, 100), Q(1, 4)), '5': (Q(2, 100), Q(3, 100), Q(1, 2))}
res = {}
for k, (fe2, fe3, ag) in cases.items():
    q = fe3/(fe2*ag)
    res[k] = 'Ag' if q < K else ('Pt' if q > K else 'none')
    print(k, float(q), res[k])
assert res == {'2': 'Ag', '3': 'Pt', '4': 'Pt', '5': 'none'}
for e in ['Fe^3+ + e- -> Fe^2+', 'Ag+ + e- -> Ag', 'Ag -> Ag+ + e-']:
    assert check_equation(e) is None
print('THE-EQL-032 OK')
