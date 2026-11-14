# THE-ELC-054 検証：金属 1.00 kg を電解で得るのに要する電気量・電力量の比較（Al・Mg・Cu）
import sympy as sp
F = 96500
data = {'Al': (sp.Rational(270, 10), 3, sp.Rational(42, 10)), 'Mg': (sp.Rational(243, 10), 2, sp.Rational(50, 10)), 'Cu': (sp.Rational(635, 10), 2, sp.Rational(25, 100))}
res = {}
for m, (M, z, V) in data.items():
    Q = 1000/M*z*F; E = Q*V; kWh = E/1000/3600
    res[m] = (Q, kWh); print(m, 'Q=%.3e C' % float(Q), 'E=%.3g kWh' % float(kWh))
assert f"{float(res['Al'][0]):.2e}" == '1.07e+07' and f"{float(res['Al'][1]):.3g}" == '12.5'
assert f"{float(res['Mg'][0]):.2e}" == '7.94e+06' and f"{float(res['Mg'][1]):.3g}" == '11'
assert f"{float(res['Cu'][0]):.2e}" == '3.04e+06' and f"{float(res['Cu'][1]):.2g}" == '0.21'
ratio = res['Al'][1]/res['Cu'][1]; qr = res['Al'][0]/res['Cu'][0]
print('Al/Cu energy', float(ratio), ' charge ratio', float(qr), ' voltage ratio', 4.2/0.25)
assert round(float(ratio)) == 59 and round(float(qr), 2) == 3.53
