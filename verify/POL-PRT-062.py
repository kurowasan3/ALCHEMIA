# POL-PRT-062 検証：Lys-Lys-Glu の等電点
from sympy import symbols, nsolve, log, Rational as Q
pH = symbols('pH')
def charge(groups):
    c = 0
    for pk, t in groups:
        c += 1/(1 + 10**(pH - pk)) if t == 'b' else -1/(1 + 10**(pk - pH))
    return c
P = [(3.2, 'a'), (4.3, 'a'), (7.9, 'b'), (10.6, 'b'), (10.6, 'b')]
D = [(3.2, 'a'), (4.3, 'a'), (7.9, 'b'), (10.6, 'b')]
pI_P = float(nsolve(charge(P), pH, 9)); pI_D = float(nsolve(charge(D), pH, 6))
approx = (7.9 + 10.6)/2 - 0.30/2
print(pI_P, approx, pI_D)
assert round(pI_P, 1) == 9.1 and round(approx, 2) == 9.10 and round(pI_D, 1) == 6.1 and round((4.3+7.9)/2, 1) == 6.1
h = 10**-approx; print('h/K1', h/10**-7.9, 'K2/h', 10**-10.6/h)
naoh = float(3 - charge(P).subs(pH, 7.9)); print('NaOH eq at 7.9:', naoh)
assert round(naoh, 1) == 2.5
assert float(charge(P).subs(pH, 0)) > 2.99 and float(charge(P).subs(pH, 14)) < -1.99
print('OK')
