# THE-ION-095 検証：塩酸＋Cu2+ の NaOH 滴定（沈殿の生成と pH）
import math
from scipy.optimize import brentq
Ksp, Kw = 2.2e-20, 1.0e-14
nH, nCu, V0, c = 2.0e-3, 1.0e-3, 0.050, 0.100
def state(Vml):
    V = Vml/1000; Vt = V0+V; nOH = c*V
    # 電荷収支を [H+] について解く（沈殿 s mol は [Cu2+][OH-]^2 = Ksp で決まる，s>=0）
    def f(h):
        oh = Kw/h; cu_sat = Ksp/oh**2
        cu = min(nCu/Vt, cu_sat)
        # 2[Cu2+] + [H+] + [Na+] = [Cl-] + [OH-]
        return 2*cu + h + nOH/Vt - (nH + 2*nCu)/Vt - oh
    h = brentq(f, 1e-14, 1.0, xtol=1e-30, rtol=1e-14)
    oh = Kw/h; cu = min(nCu/Vt, Ksp/oh**2)
    return -math.log10(h), 1 - cu*Vt/nCu
for v in (19.0, 19.9, 19.99, 20.0, 20.01, 25.0, 30.0, 35.0, 39.96, 40.0, 40.1):
    pH, frac = state(v); print(v, round(pH, 3), 'precip frac', round(frac, 5))
# 解答の値（単純化した計算）
pH_on = 14 + math.log10(math.sqrt(Ksp/(nCu/0.070)))
pH30 = 14 + math.log10(math.sqrt(Ksp/(5e-4/0.080)))
pH999 = 14 + math.log10(math.sqrt(Ksp/(1e-6/0.090)))
d = 0.5*math.log10((7.5e-4/0.075)/(2.5e-4/0.085))
print('onset', pH_on, 'half', pH30, '99.9%', pH999, 'dpH 25-75', d)
assert round(pH_on, 2) == 5.09 and round(pH30, 2) == 5.27 and round(pH999, 1) == 6.6 and round(d, 2) == 0.27
assert abs(state(30.0)[0]-pH30) < 0.01 and abs(state(25.0)[0]-pH_on) < 0.2
# 酸の当量点の 0.1 mL 前の pH < 沈殿開始 pH
assert state(19.9)[1] < 1e-9 and state(19.9)[0] < pH_on
# 当量点直前（0.01 mL 前）でも沈殿しない
assert state(19.99)[1] < 1e-9
KspM = 6.25e-3*(10**(7.53-14))**2; print('KspM', KspM); assert round(KspM*1e16, 1) == 7.2
print('OK')
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['Cu^2+ + 2OH- -> Cu(OH)2', 'H+ + OH- -> H2O', 'Cu(OH)2 <=> Cu^2+ + 2OH-']:
    assert check_equation(e) is None
