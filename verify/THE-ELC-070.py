# THE-ELC-070 鉛蓄電池の直列・並列と電解液の変化
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, N
for e in ['Pb + SO4^2- -> PbSO4 + 2e-', 'PbO2 + 4H+ + SO4^2- + 2e- -> PbSO4 + 2H2O',
          'Pb + PbO2 + 2H2SO4 -> 2PbSO4 + 2H2O', 'Cu^2+ + 2e- -> Cu', '2H2O -> O2 + 4H+ + 4e-']:
    assert check_equation(e) is None, e
Cu = R(635, 10); M_H2SO4 = 98; M_H2O = 18
nCu = R(381, 100) / Cu          # 3.81 g Cu
ne = 2 * nCu                    # 外部回路の電子
print('n(e-) total =', N(ne))
assert ne == R(12, 100)
w0, x0 = 200, R(300, 1000)      # 各電池 200 g, 30.0 %
def after(ne_cell):
    s = w0 * x0 - ne_cell * M_H2SO4
    m = w0 - ne_cell * M_H2SO4 + ne_cell * M_H2O
    return N(100 * s / m, 6)
cA = after(ne)          # 6 直列: 各電池に全電子
cB = after(ne / 2)      # 3 直列 x 2 並列: 各電池に半分
print('A:', cA, ' B:', cB)
assert round(float(cA), 1) == 25.3 and round(float(cB), 1) == 27.7
# 負極の質量増加（B の 1 個）: Pb -> PbSO4 で SO4 = 96 g/mol per 2e-
dneg = (ne / 2) / 2 * 96
print('neg electrode +', N(dneg))
assert dneg == R(288, 100)
# 全体の H2SO4 消費
totA = 6 * ne; totB = 6 * ne / 2
print('total H2SO4 A,B =', N(totA), N(totB))
assert totA == 2 * totB
print('all OK')
