# THE-THC-050 検証：リン酸と NaOH の混合熱（混合の順序によらない）
from sympy import Rational as R, nsimplify
c=4.2  # J/(g K), 密度 1.0 g/mL
# 実験1：H3PO4 0.10 mol/L 50 mL (5.0 mmol) に NaOH 0.10 mol/L を加える
n=5.0e-3
Q50 = 100*c*0.74      # J
Q100= 150*c*0.92
q1 = Q50/n/1000; q2=(Q100-Q50)/n/1000
print('Q50=%.1f J Q100=%.1f J q1=%.1f q2=%.1f kJ/mol'%(Q50,Q100,q1,q2))
assert round(q1)==62 and round(q2)==54
# 実験2(a)：NaOH 100 mL (10 mmol) に H3PO4 50 mL (5 mmol)→ 全部 HPO4^2- （実験1の100 mL 点と同じ終状態）
Qa = 5e-3*(q1+q2)*1000
print('Qa=%.1f J dT=%.2f K'%(Qa, Qa/(150*c)))
assert abs(Qa-Q100)<1e-9 and round(Qa/(150*c),2)==0.92
# 実験2(b)：NaOH 100 mL に H3PO4 100 mL (10 mmol)→ 全部 H2PO4-
Qb = 10e-3*q1*1000
print('Qb=%.1f J (%.2f kJ) dT=%.2f K'%(Qb,Qb/1000, Qb/(200*c)))
assert round(Qb/1000,2)==0.62 and round(Qb/(200*c),2)==0.74
# 実験2(c)：NaOH 100 mL に H3PO4 25 mL (2.5 mmol)→ PO4^3- 2.5 mmol + OH- 2.5 mmol 残る：第3段の熱が必要
print('(c) needs q3: PO4^3- formed', 2.5e-3)
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['H3PO4 + OH- -> H2PO4- + H2O','H2PO4- + OH- -> HPO4^2- + H2O','HPO4^2- + OH- -> PO4^3- + H2O','H3PO4 + 2OH- -> HPO4^2- + 2H2O','PO4^3- + H3PO4 -> HPO4^2- + H2PO4-', 'PO4^3- + 2H3PO4 -> 3H2PO4-']:
    print(e, check_equation(e))
