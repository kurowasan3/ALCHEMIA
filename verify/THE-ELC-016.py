import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-016/item.json'
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

ne = Q(635, 1000)/Q(635, 10)*2; close(ne, 0.0200)
nZn = Q(327, 1000)/Q(654, 10); close(nZn, 0.005)
close(2*nZn/ne*100, 50.0)
H2 = (ne - 2*nZn)/2; O2 = ne/4
close((H2+O2)*22.4e3, 224)
for e in ['Cu -> Cu^2+ + 2e-', 'Zn^2+ + 2e- -> Zn', '2H2O + 2e- -> H2 + 2OH-', '2H2O -> O2 + 4H+ + 4e-']:
    assert check_equation(e) is None
print('THE-ELC-016 OK')
