import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-PER-013/item.json'
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

O, P, S, Cl = 16, 31, 32, Q(355, 10)
oxides = {'P4O10': (4*P + 10*O, 12), 'SO3': (S + 3*O, 2), 'Cl2O7': (2*Cl + 7*O, 2)}
hit = []
for k, (M, nb) in oxides.items():
    v = Q(142, 100)/M*nb; print(k, M, float(v))
    if abs(float(v) - 0.06) < 1e-4: hit.append(k)
assert hit == ['P4O10']
for e in ['P4O10 + 6H2O -> 4H3PO4', 'SO3 + H2O -> H2SO4', 'Cl2O7 + H2O -> 2HClO4']:
    assert check_equation(e) is None
print('THE-PER-013 OK')
