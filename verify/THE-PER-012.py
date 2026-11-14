import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/THE-PER-012/item.json'
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

el = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json'))
el = el if isinstance(el, list) else el['elements']
E = {e['sym']: e for e in el}
assert max(el, key=lambda e: e.get('en') or 0)['sym'] == 'F'
p3 = [e for e in el if e['period'] == 3]
assert max(p3, key=lambda e: e['ie1'])['sym'] == 'Ar'
cand = ['F', 'Mg', 'N', 'Si', 'P', 'Cl', 'S']
assert max(cand, key=lambda s: E[s]['radius']) == 'Mg'
print('THE-PER-012 OK')
