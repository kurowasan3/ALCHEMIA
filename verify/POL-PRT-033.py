# POL-PRT-033 検証：アミノ酸のエチルエステル塩酸塩の質量からアミノ酸を決める
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, nsimplify
M = {'H': R(1), 'C': R(12), 'N': R(14), 'O': R(16), 'Cl': R(71, 2)}
def mass(f):
    import re
    return sum(M[e]*(int(n) if n else 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', f))
aa = {'Gly': 'C2H5NO2', 'Ala': 'C3H7NO2', 'Ser': 'C3H7NO3', 'Asp': 'C4H7NO4', 'Glu': 'C5H9NO4', 'Lys': 'C6H14N2O2', 'Phe': 'C9H11NO2'}
Mw = {k: mass(v) for k, v in aa.items()}
assert (Mw['Ala'], Mw['Glu'], Mw['Lys']) == (89, 147, 146)
# 生成物の式（エステル化 a 個，塩酸塩 b 個）と反応式
prod = {'Ala': ('C5H12ClNO2', 1, 1), 'Glu': ('C9H18ClNO4', 2, 1), 'Lys': ('C8H20Cl2N2O2', 1, 2)}
for k, (f, a, b) in prod.items():
    lhs = f"{aa[k]} + {a if a > 1 else ''}C2H5OH + {b if b > 1 else ''}HCl"
    rhs = f"{f} + {a if a > 1 else ''}H2O"
    e = f"{lhs} -> {rhs}"
    assert check_equation(e) is None, e
    assert mass(f) == Mw[k] + 28*a + R(73, 2)*b
# (1) アラニン 13.35 g
n = R(1335, 100)/89
p1 = n*mass('C5H12ClNO2'); print('(1)', n, float(p1)); assert abs(p1 - R(23025, 1000)) < R(1, 10**6)
# (2) X 14.6 g → 24.7 g
for k in ('Lys', 'Glu'):
    a, b = prod[k][1], prod[k][2]
    print('(2)', k, float(R(146, 10)/Mw[k]*mass(prod[k][0])))
assert round(float(R(146, 10)/146*mass('C8H20Cl2N2O2')), 1) == 24.7
assert round(float(R(146, 10)/147*mass('C9H18ClNO4')), 1) == 23.8
# (3) Y 0.0400 mol，5.88 g，質量増加 3.70 g → 1 mol あたり 92.5 = 28a + 36.5b の整数解
g = R(370, 100)/R(4, 100)
sols = [(a, b) for a in range(0, 6) for b in range(0, 6) if 28*a + R(73, 2)*b == g]
print('(3)', g, sols, R(588, 100)/R(4, 100)); assert sols == [(2, 1)] and R(588, 100)/R(4, 100) == 147
print('POL-PRT-033 OK')
