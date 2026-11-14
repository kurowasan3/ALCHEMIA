import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-ELC-019/item.json'
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

ne = 2*Q(3175, 1000)/Q(635, 10); close(ne, 0.1)
cu_diss = ne/2 - Q(5, 1000)
nFe = (ne - 2*cu_diss)/2; mFe = nFe*56; close(mFe, 0.28)
mAg = Q(330, 100) - cu_diss*Q(635, 10) - mFe
print('Ag', float(mAg)); assert sig(mAg, 2) == 0.16
for e in ['Cu^2+ + 2e- -> Cu', 'Fe -> Fe^2+ + 2e-', 'Pb^2+ + SO4^2- -> PbSO4']:
    assert check_equation(e) is None
print('THE-ELC-019 OK')
