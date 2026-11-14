import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, N
# 25 mL aliquot
n_NaOH = R(100,1000)*R(400,10000)          # 0.100 * 40.0 mL
n_CO3 = R(100,1000)*R(60,10000)/2          # (46.0-40.0) mL / 2
w_CO3 = n_CO3*10*106/2*100; w_NaOH = n_NaOH*10*40/2*100
print('Na2CO3 %', N(w_CO3,4), 'NaOH %', N(w_NaOH,4)); assert round(float(w_CO3),1)==15.9 and round(float(w_NaOH),1)==80.0
assert float(n_CO3*10*106 + n_NaOH*10*40) < 2.00   # NaCl 分が正
# MO titre if BaCO3 left: equals total base -> 46.0 mL
assert (n_NaOH + 2*n_CO3)/R(100,1000)*1000 == 46
# (4) 2NaOH + CO2 -> Na2CO3: total HCl equivalent unchanged
x = R(1,1000); assert (n_NaOH-2*x) + 2*(n_CO3+x) == n_NaOH + 2*n_CO3
assert not check_item('THE-ACB-014')
