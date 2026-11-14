import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-024/item.json'
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

ne = Q(9650)/96500; close(ne, 0.1)
t = symbols('t')
# cathode HCl change = -ne*(1-t)
assert (-ne*(1 - 1)) == 0
after = Q(1)*Q(241, 10)/Q(25)*Q(1, 2); dec = Q(1, 2) - after; close(dec, 0.018)
tt = solve(ne*(1 - t) - dec, t)[0]; close(tt, 0.82)
close(ne*tt, 0.082)
r = tt/(1 - tt); print('ratio', float(r)); assert sig(r, 2) == 4.6
ca = (Q(1, 2) - ne*tt)/Q(1, 2); print('anode', float(ca)); assert sig(ca, 2) == 0.84
for e in ['2H+ + 2e- -> H2', '2Cl- -> Cl2 + 2e-']:
    assert check_equation(e) is None
print('THE-ELC-024 OK')
