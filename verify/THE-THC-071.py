# THE-THC-071 C1 化合物の燃焼エンタルピーと酸化の段階（ヘスの法則）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
# 燃焼エンタルピー（25 ℃，生成物 CO2(気)・H2O(液)），kJ/mol：文献値 CH4 -890.8, CH3OH(l) -726.1, HCHO(g) -570.7, HCOOH(l) -254.6, CO -283.0
c = {'CH4': -891, 'CH3OH': -726, 'HCHO': -571, 'HCOOH': -255, 'CO': -283}
steps = [('2CH4 + O2 -> 2CH3OH', 'CH4', 'CH3OH'), ('2CH3OH + O2 -> 2HCHO + 2H2O', 'CH3OH', 'HCHO'),
         ('2HCHO + O2 -> 2HCOOH', 'HCHO', 'HCOOH')]
for eq, a, b in steps:
    print(eq, check_equation(eq), ' ΔH per mol =', c[a] - c[b])
print('HCOOH -> CO + H2O', check_equation('HCOOH -> CO + H2O'), ' ΔH =', c['HCOOH'] - c['CO'])
for eq in ['CH4 + 2O2 -> CO2 + 2H2O', '2CH3OH + 3O2 -> 2CO2 + 4H2O', 'HCHO + O2 -> CO2 + H2O', '2HCOOH + O2 -> 2CO2 + 2H2O', '2CO + O2 -> 2CO2']:
    print(eq, check_equation(eq))
# 炭素の酸化数
print('ox numbers: CH4 -4, CH3OH -2, HCHO 0, HCOOH +2, CO +2, CO2 +4')
# 和の確認：CH4 の燃焼 = 3段階 + HCOOH の燃焼
print('sum check', (c['CH4']-c['CH3OH'])+(c['CH3OH']-c['HCHO'])+(c['HCHO']-c['HCOOH'])+c['HCOOH'], '=', c['CH4'])
