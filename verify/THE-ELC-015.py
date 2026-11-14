import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-015/item.json'
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

F = 96500
Qtot = Q(1, 2) * (32*60 + 10); close(Qtot, 965)
nI = Q(648, 1000) / 108; QI = nI * F; close(QI, 579)
QII = Qtot - QI; close(QII, 386)
ne = QII / F
gas = ne/2 + ne/4; close(gas*22.4e3, 67.2)
# parallel circuit
V, R, RI, RII = symbols('V R R_I R_II', positive=True)
Rp = RI*RII/(RI+RII); It = V/(R+Rp); II = It*Rp/RI
base = dict(V=10, R=5, RI=4, RII=6)
f = lambda **k: float(II.subs({V: k['V'], R: k['R'], RI: k['RI'], RII: k['RII']}))
b0 = f(**base)
assert f(**{**base, 'RII': 3}) < b0      # (4) conc. up -> R_II down -> decrease
assert f(**{**base, 'RI': 2}) > b0       # (5) closer -> R_I down -> increase
assert f(**{**base, 'RII': 1e12}) > b0   # (6) remove plate -> increase
assert f(**{**base, 'R': 8}) < b0        # (7) R up -> decrease
for e in ['Ag+ + e- -> Ag', '2H2O + 2e- -> H2 + 2OH-', '2H2O -> O2 + 4H+ + 4e-']:
    assert check_equation(e) is None
print('THE-ELC-015 OK')
