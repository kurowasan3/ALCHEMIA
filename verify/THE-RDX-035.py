# THE-RDX-035 検証：二酸化硫黄水溶液のヨウ素による酸化と終点の pH
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as Q, log, N
assert check_equation('SO2 + I2 + 2H2O -> H2SO4 + 2HI') is None
assert check_equation('SO2 + H2O2 -> H2SO4') is None
c, V0, H = Q(5, 10000), Q(50), Q(1, 1000)    # mol/L, mL, mol/L (pH 3.0)
y = symbols('y', positive=True)
ys = solve(4*c*y - H*(V0 + y), y)[0]; print('(3) y =', ys); assert ys == 50
nSO2 = c*ys/1000; conc = nSO2/(V0/1000); print('   n(SO2) =', nSO2, ' conc =', conc)
assert nSO2 == Q(25, 1000000) and conc == Q(5, 10000)
# HSO4- の電離（Ka2 = 1.0e-2）を考えたときの補正の大きさ（解説用）
Ka2 = 1.0e-2; Hf = 1.0e-3
frac_SO4 = Ka2/(Ka2 + Hf); print('   SO4^2- fraction at pH 3 =', frac_SO4, ' H+ per S =', 1 + frac_SO4 + 2)
# (4) H2O2 で酸化した場合：H+ は S あたり 2（体積は 50 mL のまま、加えた H2O2 の体積は無視）
H2 = 2*nSO2/(V0/1000); print('(4) [H+] =', H2, ' pH =', N(-log(H2, 10), 3)); assert H2 == Q(1, 1000)
# 同じ試料をヨウ素で酸化したあと（全量 100 mL）は [H+] = 1.0e-3 で同じ pH 3.0 → 一致する
print('OK')
