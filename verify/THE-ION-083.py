# THE-ION-083 検証：霧粒・雨粒への二酸化硫黄の取り込みと硫酸水素イオンの電離平衡
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, sqrt, log, N, Rational as Q
for e in ['SO2 + H2O <=> H2SO3', 'H2SO3 <=> H+ + HSO3-', '2HSO3- + O2 -> 2HSO4-',
          'H2SO4 -> H+ + HSO4-', 'HSO4- <=> H+ + SO4^2-']:
    assert check_equation(e) is None, e
n_air = 44.6; n_S = 4.0e-7; P = 1.0e5            # 1 m3 の空気（mol），SO2（mol），全圧 Pa
kH = 1.2e-5      # mol/(L Pa)  H2SO3（溶けた SO2）の溶解度
Ka1 = 1.7e-2; Ka2 = 1.0e-2
# (1) 酸化が起こらない場合、pH 5.0 の霧粒（1 m3 中 0.10 g = 1.0e-4 L）に入る割合
L = 1.0e-4
pSO2 = n_S/n_air*P                                   # 気相にほぼ全量が残るとした分圧
c_tot = kH*pSO2*(1 + Ka1/1.0e-5)                     # H2SO3 + HSO3-
frac = c_tot*L/n_S
print('(1) p(SO2) =', pSO2, 'Pa  dissolved conc =', c_tot, ' fraction =', frac)
assert f'{frac:.1e}' == '4.6e-03'
# 自己無撞着な pH（[H+] = [HSO3-] として）
Hsc = sqrt(Ka1*kH*pSO2); print('   self-consistent [H+] =', N(Hsc, 3), ' pH =', N(-log(Hsc, 10), 3), ' fraction =', N(Hsc*L/n_S, 2))
assert float(Hsc*L/n_S) < 0.01
# (2) 酸化が進むとき：K3 = [HSO4-]^2/(p(O2)[HSO3-]^2) = 1e80 (atm^-1 相当)
ratio = 1/sqrt(1e80*0.21); print('(2) [HSO3-]/[HSO4-] =', N(ratio, 3))
# (3) 霧（L = 1.0e-4 L/m3）：C = 4.0e-3 mol/L
def H_of(C):
    x = symbols('x', positive=True)
    xs = [s for s in solve(x*(C + x) - Ka2*(C - x), x) if 0 < s < C][0]
    return C + xs
C_fog = n_S/1.0e-4; Hf = H_of(C_fog); pHf = -log(Hf, 10)
print('(3) fog C =', C_fog, ' [H+] =', N(Hf, 4), ' pH =', N(pHf, 4)); assert round(float(pHf), 1) == 2.2
# (4) 雨雲（L = 1.0e-3 L/m3）：C = 4.0e-4 mol/L
C_rain = n_S/1.0e-3; Hr = H_of(C_rain); pHr = -log(Hr, 10)
print('(4) rain C =', C_rain, ' [H+] =', N(Hr, 4), ' pH =', N(pHr, 4)); assert round(float(pHr), 1) == 3.1
print('   [H+] ratio fog/rain =', N(Hf/Hr, 3), '  2-valent approximation pH:', N(-log(2*C_fog, 10), 3), N(-log(2*C_rain, 10), 3))
print('   SO4^2- fraction fog', N((Hf - C_fog)/C_fog, 3), ' rain', N((Hr - C_rain)/C_rain, 3))
assert round(float(Hf/Hr), 1) == 8.3
# 解答で使う 2 桁の値
print('rounded:', f'{float(Hf):.1e}', f'{float(Hr):.1e}')
print('OK')
