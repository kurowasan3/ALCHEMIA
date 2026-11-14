# THE-ION-100 検証：ZnS の沈殿とヒドロキシド錯イオンの競合
import math
K1, K2, Kw, b, Ksp, C = 1e-7, 1e-14, 1e-14, 1e15, 1e-24, 1e-3
def ZnT(pH):
    h = 10**-pH; s = C*K1*K2/(h*h + K1*h + K1*K2); z = Ksp/s
    return s, z*(1 + b*(Kw/h)**4)
s3, z3 = ZnT(3.0); print(s3, z3); assert abs(s3/1e-18 - 1) < 1e-3 and abs(z3/1e-6 - 1) < 1e-3
z9, z12 = ZnT(9.0)[1], ZnT(12.0)[1]; print(z9, z12)
assert round(z9/1e-16, 1) == 1.0 and round(z12/1e-12, 1) == 1.0
best = min((ZnT(6 + i*1e-4)[1], 6 + i*1e-4) for i in range(80001)); print(best)
approx = (41 - 0.48)/4; print(approx); assert abs(best[1] - approx) < 0.01 and round(approx, 1) == 10.1
h = 10**-approx; assert h < K1/100 and h > K2*1e3
print('Zn(OH)2 IP', (Ksp/(C*K2/h))*(Kw/h)**2)
print('OK')
