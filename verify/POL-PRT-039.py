# POL-PRT-039 検証：修飾したアミノ酸の pH 6.0・pH 1.0 での正味の電荷（電気泳動の向き）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
# 代表的な文献値の pKa（酸性基: COOH，塩基性基: NH3+ の pKa）
DATA = {
 'グリシン':               {'acid': [2.34], 'base': [9.60]},
 'N-アセチルグリシン':      {'acid': [3.67], 'base': []},
 'グリシンエチルエステル':   {'acid': [], 'base': [7.7]},
 'リシン':                 {'acid': [2.18], 'base': [8.95, 10.53]},
 'グルタミン酸':            {'acid': [2.19, 4.25], 'base': [9.67]},
 'Nα-アセチルリシン':       {'acid': [3.5], 'base': [10.5]},
 'Nα,Nε-ジアセチルリシン':  {'acid': [3.5], 'base': []},
}
def charge(d, pH):
    q = 0.0
    for pk in d['acid']: q -= 1 / (1 + 10 ** (pk - pH))
    for pk in d['base']: q += 1 / (1 + 10 ** (pH - pk))
    return q
def direction(q):
    return '陰極' if q > 0.3 else '陽極' if q < -0.3 else 'ほとんど動かない'
expect6 = {'グリシン': 'ほとんど動かない', 'N-アセチルグリシン': '陽極', 'グリシンエチルエステル': '陰極',
           'リシン': '陰極', 'グルタミン酸': '陽極', 'Nα-アセチルリシン': 'ほとんど動かない', 'Nα,Nε-ジアセチルリシン': '陽極'}
expect1 = {'グリシン': '陰極', 'N-アセチルグリシン': 'ほとんど動かない', 'グリシンエチルエステル': '陰極',
           'リシン': '陰極', 'グルタミン酸': '陰極'}
for nm, d in DATA.items():
    q6, q1 = charge(d, 6.0), charge(d, 1.0)
    print(f'{nm}: pH6 {q6:+.2f} {direction(q6)}   pH1 {q1:+.2f} {direction(q1)}')
    assert direction(q6) == expect6[nm], nm
    if nm in expect1: assert direction(q1) == expect1[nm], nm
for e in ['H2NCH2COOH + (CH3CO)2O -> CH3CONHCH2COOH + CH3COOH',
          'H2NCH2COOH + C2H5OH -> H2NCH2COOC2H5 + H2O']:
    r = check_equation(e); print(e, r); assert r is None
print('ALL OK')
