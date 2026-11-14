import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/INO-GPR-029/item.json'
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

import itertools
g = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json'))['gases']
gases = ['H2', 'O2', 'CO', 'CO2', 'NO', 'NO2', 'NH3', 'HCl', 'H2S', 'Cl2']
col = {k: g[k]['color'] != '無色' for k in gases}
odor = {k: g[k]['smell'] != '無臭' for k in gases}
combust = {'H2': {'H2O'}, 'CO': {'CO2'}, 'H2S': {'H2O', 'SO2'}, 'NH3': {'N2', 'H2O'}}
sols = []
for perm in itertools.permutations(gases, 8):
    A, B_, C, D, E_, F, G, H = perm
    if any(col[x] for x in perm[:7]) or not col[H]: continue
    if any(odor[x] for x in perm[:5]) or not all(odor[x] for x in (F, G, H)): continue
    if E_ != 'NO': continue
    if combust.get(A) != {'H2O'} or combust.get(C) != {D}: continue
    if D != 'CO2' or B_ != 'O2': continue
    if g[F]['nature'] != '塩基性': continue
    if G != 'H2S': continue
    if H != 'Cl2': continue   # bleaching aqueous solution (HClO) only for Cl2 among colored
    sols.append(perm)
print(sols); assert sols == [('H2', 'O2', 'CO', 'CO2', 'NO', 'NH3', 'H2S', 'Cl2')]
assert '赤褐' in g['NO2']['color'] and '酢酸鉛' in g['H2S']['test']
for e in ['2CO + O2 -> 2CO2', '2NO + O2 -> 2NO2', 'Cl2 + H2O -> HCl + HClO']:
    assert check_equation(e) is None
print('INO-GPR-029 OK')
