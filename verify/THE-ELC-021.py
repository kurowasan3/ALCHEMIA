import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-021/item.json'
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

# simulate compartments A..E with membranes X(cation) Y(anion) X Y between them
mem = ['X', 'Y', 'X', 'Y']
na = [0]*5; cl = [0]*5; oh = [0]*5
# per 1 mol e-: anode in 0 oxidizes Cl-, cathode in 4 produces OH-
cl[0] -= 1; oh[4] += 1
for i in range(5):
    # cation tries to move right (to i+1), anion left (to i-1); 1 mol of charge must cross each membrane
    pass
for j, m in enumerate(mem):  # membrane between j and j+1
    if m == 'X':
        na[j] -= 1; na[j+1] += 1
    else:
        cl[j+1] -= 1; cl[j] += 1
nacl = [min(0, 0) for _ in range(5)]
# NaCl change = change in Na+ paired with Cl- (E: Na+ unchanged, Cl- -1 replaced by OH-)
for i in range(5):
    charge = na[i] - cl[i] - oh[i]
    assert charge == 0, (i, charge)
delta = [cl[i] for i in range(5)]
print('dNaCl per e-', delta)
assert [i for i in range(5) if delta[i] > 0] == [1, 3]
nH2 = Q(100000)*Q(249, 100)/(Q(8310)*300); ne = 2*nH2; print('ne', float(ne)); assert sig(ne, 2) == 0.20
ne = Q(2, 10)
I = ne*96500/7200; print('I', float(I)); assert sig(I, 2) == 2.7
pH = 14 - (1 - 0.30); assert abs(pH - 13.3) < 1e-9
conc = [Q(1, 2) + delta[i]*ne for i in range(5)]
print([float(c) for c in conc]); assert [c*10 for c in conc] == [3, 7, 3, 7, 3]
for e in ['2Cl- -> Cl2 + 2e-', '2H2O + 2e- -> H2 + 2OH-']:
    assert check_equation(e) is None
print('THE-ELC-021 OK')
