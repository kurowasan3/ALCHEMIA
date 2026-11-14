# THE-BND-061 検証：ホタル石型酸化物（CeO2）への異価イオン添加と酸化物イオン空孔・酸素ポンプ
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
R = sp.Rational
# 単位格子：陽イオン 8*1/8+6*1/2=4，陰イオン 8
cat = 8*R(1,8)+6*R(1,2); an = 8; assert (cat, an) == (4, 8)
def vac(dopant_charge, x):  # 陽イオン1個あたりの O の数と空孔の割合
    O = ((1-x)*4 + x*dopant_charge)/R(2)
    return O, (2-O)/2
# CeO2 : Gd2O3 = 0.80 : 0.10 → 陽イオン Ce 0.80, Gd 0.20
x = R(20,100)
O, f = vac(3, x); print('Gd: O per cation', O, 'vacancy fraction', f); assert f == R(5,100)
print('vacancies per cell', f*an); assert f*an == R(4,10)
O2, f2 = vac(2, x); print('Ca: vacancy fraction', f2); assert f2 == R(10,100)
# 酸素ポンプ
q = R(50,100)*386; ne = q/96500; nO2 = ne/4
V = nO2*R(224,10)*1000
print('e-', ne, 'O2', nO2, 'V mL', float(V)); assert round(float(V)) == 11
for e in ['O2 + 4e- -> 2O^2-', '2O^2- -> O2 + 4e-']:
    assert check_equation(e) is None
print('OK')
