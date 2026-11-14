# THE-THC-065 オゾンの熱化学 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N
BO2=498; AO3=604   # O2 の結合エネルギー, O3 の原子化エネルギー (kJ/mol)
dH1 = 3*BO2 - 2*AO3          # 3O2 -> 2O3
dH2 = BO2 - AO3              # O + O2 -> O3
dH3 = AO3 - 2*BO2            # O + O3 -> 2O2
dH4 = AO3 - BO2              # O3 -> O2 + O (必要なエネルギー)
avg = AO3/2                  # O3 中の O-O 結合 1 本あたり
print(dH1, dH2, dH3, dH4, avg)
assert (dH1,dH2,dH3,dH4)==(286,-106,-392,106) and avg==302
# 大気中のオゾン 1.00 m2 あたり：0 ℃, 1.013e5 Pa で厚さ 3.00 mm
P=1.013e5; V=3.00; R=8.31e3; T=273
n=P*V/(R*T); print('n O3 =', n)
assert abs(n-0.1340)<5e-4
q = n*dH1/2      # 2O3 -> 3O2 で放出 (kJ) ; O3 1 mol あたり 143 kJ
print('heat =', q); assert round(q,1)==19.2
for e in ['3O2 -> 2O3','O + O2 -> O3','O + O3 -> 2O2','O3 -> O2 + O','O2 -> 2O','2O3 -> 3O2']:
    assert check_equation(e) is None, e
print('OK')
