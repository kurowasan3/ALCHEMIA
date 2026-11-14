# ORG-ALC-031 検証：カルボニル化合物の水和平衡・遊離カルボニルの割合・18O 交換
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
K = {'クロラール': sp.Rational(28000), 'プロパナール': sp.Rational(85, 100), '2-メチルプロパナール': sp.Rational(6, 10), 'アセトン': sp.Rational(14, 10000)}
# (1) 電子求引性（Cl3C-）で増加，アルキル基の数・かさ高さで減少：順序の確認
order = sorted(K, key=lambda k: -K[k]); print(order)
assert order == ['クロラール', 'プロパナール', '2-メチルプロパナール', 'アセトン']
# (2) 遊離カルボニルの割合 f = 1/(1+K)
f = {k: 1/(1+v) for k, v in K.items()}
print('f(プロパナール)=', float(f['プロパナール']), ' f(クロラール)=', float(f['クロラール']))
assert round(float(f['プロパナール']), 2) == 0.54
assert f"{float(f['クロラール']):.1e}" == '3.6e-05'
# (3) 初速度の比 v(クロラール)/v(プロパナール) = 100 × f比
r = 100 * f['クロラール'] / f['プロパナール']
print('rate ratio', float(r), ' inverse', float(1/r)); assert f"{float(r):.1e}" == '6.6e-03'
# (4) 18O の分配：O の総量 1.0 + 9.0 mol，18O は 0.90 mol → 統計的に均一
x = sp.Rational(90, 100) / (1 + 9)
print('18O fraction', float(x)); assert float(x) == 0.09
for e in ['CH3CH2CHO + H2O <=> CH3CH2CH(OH)2', 'CCl3CHO + H2O <=> CCl3CH(OH)2', 'CH3COCH3 + H2O <=> (CH3)2C(OH)2']:
    assert check_equation(e) is None, e
print('equations OK')
