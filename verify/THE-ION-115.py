# THE-ION-115 検証：リン酸カルシウムの溶解度の pH 依存性と臨界 pH
import sympy as sp
K1, K2, K3 = sp.Float('7.5e-3'), sp.Float('6.2e-8'), sp.Float('4.8e-13')
Ksp = sp.Float('2.0e-29')
def a3(h): return K1*K2*K3/(h**3 + K1*h**2 + K1*K2*h + K1*K2*K3)
res = {}
for pH in (7, 5):
    h = sp.Float(10)**(-pH); a = a3(h); S = (Ksp/(108*a**2))**sp.Rational(1, 5)
    res[pH] = (a, S); print('pH', pH, 'alpha3 = %.3g' % a, 'S = %.3g' % S)
assert '%.2g' % res[7][0] == '1.8e-06' and '%.2g' % res[5][0] == '3e-10'
assert '%.2g' % res[7][1] == '0.00014' and '%.2g' % res[5][1] == '0.0046'
print('ratio S5/S7 = %.3g' % (res[5][1]/res[7][1]))
# K1 >> h >> K2 の領域：alpha3 ≈ K2K3/h^2，S ∝ h^(4/5)
print('per pH unit factor = %.3g' % (10**0.8))
# 臨界 pH
Ca, CP = sp.Float('1.5e-3'), sp.Float('4.0e-3')
h = sp.symbols('h', positive=True)
x = sp.symbols('x')
pHc = sp.nsolve(sp.log(Ca**3*(CP*a3(sp.Float(10)**(-x)))**2/Ksp, 10), x, (5, 7), solver='bisect'); print('critical pH (exact) = %.4f' % pHc)
# 近似（h << K1）：h(h+K2) = K2K3/alpha
alpha_need = sp.sqrt(Ksp/(Ca**3*CP**2)); print('alpha needed = %.4g' % alpha_need)
hh = sp.nsolve(h**2 + K2*h - K2*K3/alpha_need, h, 1e-6); print('approx pH = %.4f' % (-sp.log(hh, 10)))
hh2 = sp.sqrt(K2*K3/alpha_need); print('cruder (h>>K2) pH = %.4f' % (-sp.log(hh2, 10)))
assert '%.2g' % pHc == '5.9'
for pH in (6.5, 5.5):
    hv = sp.Float(10)**(-sp.Float(pH)); print(pH, 'Q/Ksp = %.3g' % (Ca**3*(CP*a3(hv))**2/Ksp))
print('OK')
