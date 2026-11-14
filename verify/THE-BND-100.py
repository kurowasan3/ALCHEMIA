# THE-BND-100 verification: Bragg condition, Miller indices for MgO (rock salt)
from sympy import sqrt, Rational, N
a = 421.2; lam = 154.0
d111 = a/float(sqrt(3)); d220 = a/float(sqrt(8))
print('d111', d111, a/1.732, 'd220', d220, a/(2*1.414))
assert round(d111) == 243 and round(a/1.732) == 243
assert round(d220) == 149 and round(a/(2*1.414)) == 149
for s_obs, s_exp, name in [(0.3656, 4, '200'), (0.5171, 8, '220')]:
    d = lam/(2*s_obs); s = (a/d)**2
    print(name, round(d, 1), round(s, 3))
    assert abs(s - s_exp) < 0.01
    assert all(abs(s - c) > 0.5 for c in (3, 4, 8, 11) if c != s_exp)
assert round(lam/(2*0.3656), 1) == 210.6 and round(lam/(2*0.5171), 1) == 148.9
assert round(a/2) == 211                       # Mg-O shortest = a/2 = d200
s2 = lam/(a/1.732); print('n=2 (111) sin', s2)
assert round(s2, 3) == 0.633
assert abs(a/1.732/2 - a/float(sqrt(12))) < 0.1  # d111/2 = d222
smax = (2*a/lam)**2; smax_mo = (2*a/71)**2; print('smax', smax, smax_mo)
assert int(smax) == 29 and 27 <= smax < 32 and smax_mo > 32
print('sin(333)', lam*float(sqrt(27))/(2*a))
assert round(lam*float(sqrt(27))/(2*a), 3) == 0.950
print('OK')
