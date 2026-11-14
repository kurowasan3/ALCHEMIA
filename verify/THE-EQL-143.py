# THE-EQL-143 Chapman 機構の定常状態 検証（数値積分で定常状態に達することも確かめる）
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
O, O3 = sp.symbols('O O3', positive=True)
j1, j3, k2, k4, O2, M = 1.0e-11, 1.0e-3, 6.0e-34, 1.0e-15, 2.0e17, 1.0e18
# 厳密な定常解
# 和をとると k4[O][O3]=j1[O2]．これを代入すると [O] の二次方程式 k2 M k4 O^2 - j1 k4 O - j1 j3 = 0
a,b,c=k2*M*k4,-j1*k4,-j1*j3
Oex=(-b+(b*b-4*a*c)**0.5)/(2*a); O3ex=j1*O2/(k4*Oex)
# 近似式
O3ap=O2*((j1*k2*M)/(j3*k4))**0.5
Oap=j3*O3ap/(k2*O2*M)
print('exact',Oex,O3ex,' approx',Oap,O3ap)
assert abs(O3ap-1.55e13)/1.55e13<0.01 and abs(Oap-1.29e8)/1.29e8<0.01
assert abs(O3ex-O3ap)/O3ap<1e-3
r=k4*Oap/j3; print('v4/v3',r); assert abs(r-1.29e-4)<0.01e-4
# 2v1 = 2v4 の確認
assert abs(j1*O2 - k4*Oex*O3ex)/(j1*O2) < 1e-6
# O2, M を 1/2
f=(0.5)*(0.5)**0.5; print('factor',f); assert round(f,3)==0.354
# 実効 k4 を 4 倍
print(1/4**0.5)
for e in ['O2 -> 2O','O3 -> O + O2','O + O3 -> 2O2','Cl + O3 -> ClO + O2','ClO + O -> Cl + O2']:
    assert check_equation(e) is None
print('OK')
