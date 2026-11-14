import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/INO-AMP-011/item.json'
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

tot = Q(448, 1)/Q(224, 10); CO = tot/4; CO2 = tot*3/4
ne = CO*2 + CO2*4; close(ne, 70)
mAl = ne/3*27; close(mAl, 630)
close(tot*12, 240)
eff = ne/(Q(100)*20*3600/96500)*100; print('eff', float(eff)); assert sig(eff, 2) == 94
close(ne/4*12, 210)
for e in ['C + O^2- -> CO + 2e-', 'C + 2O^2- -> CO2 + 4e-', 'Al^3+ + 3e- -> Al', '2H2O + 2e- -> H2 + 2OH-']:
    assert check_equation(e) is None
print('INO-AMP-011 OK')
