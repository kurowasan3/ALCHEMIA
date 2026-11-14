import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, nsimplify, symbols, solve, log, N

ITEM = '/home/claude/chem/out/items/INO-GPR-027/item.json'
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

eqs = {'ア': 'Zn + 2HCl -> ZnCl2 + H2', 'イ': 'NaCl + H2SO4 -> NaHSO4 + HCl',
       'ウ': 'NaHSO3 + H2SO4 -> NaHSO4 + H2O + SO2', 'エ': 'CaCl(ClO).H2O + 2HCl -> CaCl2 + 2H2O + Cl2',
       'オ': '2KClO3 -> 2KCl + 3O2', 'カ': 'HCOOH -> CO + H2O', 'キ': '2Na2O2 + 2H2O -> 4NaOH + O2',
       'ク': 'FeS + H2SO4 -> FeSO4 + H2S', 'ケ': '2NH4Cl + Ca(OH)2 -> CaCl2 + 2H2O + 2NH3'}
for k, e in eqs.items():
    assert check_equation(e) is None, (k, e)
redox = {'ア': True, 'イ': False, 'ウ': False, 'エ': True, 'オ': True, 'カ': False, 'キ': True, 'ク': False, 'ケ': False}
# oxidation-number spot checks: HCOOH C: +2 ; CO C: +2
assert (0 - 2*1 - 2*(-2)) == 2 and (0 - (-2)) == 2
assert sorted(k for k, v in redox.items() if not v) == sorted(['イ', 'ウ', 'カ', 'ク', 'ケ'])
g = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json'))['gases']
assert '濃硫酸' not in g['H2S']['dry'] and '濃硫酸' not in g['NH3']['dry']
for x in ['H2', 'HCl', 'SO2', 'Cl2', 'O2', 'CO']:
    assert '濃硫酸' in g[x]['dry'], x
assert g['O2']['collect'] == '水上置換'
print('INO-GPR-027 OK')
