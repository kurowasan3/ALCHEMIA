import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/INO-HAL-011/item.json'
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

x, y = symbols('x y')
s = solve([158*x + 87*y - Q(561, 100), Q(5, 2)*x + y - Q(1904, 1000)/Q(224, 10)], [x, y])
close(s[x], 0.03); close(s[y], 0.01)
p = s[x]*158/Q(561, 100)*100; print('purity', float(p)); assert sig(p, 2) == 84
for e in ['2MnO4- + 16H+ + 10Cl- -> 2Mn^2+ + 5Cl2 + 8H2O']:
    assert check_equation(e) is None
g = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json'))['gases']
assert set(g['Cl2']['dry']) == {'濃硫酸', '塩化カルシウム'}
print('INO-HAL-011 OK')
