# INO-PPT-055 検証：鉄の重量分析（Fe(OH)3 -> Fe2O3）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs = ['Fe + 2H+ -> Fe^2+ + H2',
       '2Fe^2+ + H2O2 + 2H+ -> 2Fe^3+ + 2H2O',
       'Fe^3+ + 3NH3 + 3H2O -> Fe(OH)3 + 3NH4+',
       '2Fe(OH)3 -> Fe2O3 + 3H2O',
       'NH4NO3 -> N2O + 2H2O']
for e in eqs:
    r = check_equation(e); print(e, '->', r); assert r is None
Fe, O = 55.8, 16.0
M = 2*Fe + 3*O
mFe = 0.732*2*Fe/M
pct = mFe/0.520*100
print('Fe2O3', M, 'Fe g', mFe, '%', pct)
assert f'{pct:.3g}' == '98.4'
# 強熱が不十分で FeO(OH) のまま量った場合（Fe 1 mol あたり 88.8 g vs 79.8 g）
print('FeO(OH)/ (1/2 Fe2O3)', (Fe+2*O+1.0)/(M/2))
