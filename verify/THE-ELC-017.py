import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-017/item.json'
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

for e in ['NiO(OH) + H2O + e- -> Ni(OH)2 + OH-', 'Cd + 2NiO(OH) + 2H2O -> Cd(OH)2 + 2Ni(OH)2', 'MH + NiO(OH) -> M + Ni(OH)2']:
    pass
assert check_equation('NiO(OH) + H2O + e- -> Ni(OH)2 + OH-') is None
assert check_equation('Cd + 2NiO(OH) + 2H2O -> Cd(OH)2 + 2Ni(OH)2') is None
OH = 17
close(Q(1, 10)/2*2*OH, 1.7)       # negative: +2OH per 2e-
close(Q(1, 10)*1.0, 0.10)         # positive: +H per e-
M = 59 + 16 + 17; close(M, 92)
n = Q(276, 100)/M; t = n*96500/Q(4, 10)/3600
print('hours', float(t)); assert sig(t, 2) == 2.0
print('THE-ELC-017 OK')
