# THE-ELC-085 亜鉛の湿式製錬：電解で再生する酸と硫酸イオンの収支
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R
# (1)
O2=500; H2=50
e=4*O2; z=(e-2*H2)/2; print('Zn',z,'eff',2*z/e)
# 正味の H+ = 陽極 4*O2 - 陰極 2*H2
print('net H+',4*O2-2*H2,'= 2z',2*z)
# (2) 効率80%, z=950
e2=2*950/R(8,10); print('e',e2,'O2',float(e2/4),'H2',float((e2-1900)/2),'acid',(e2-2*(e2-1900)/2)/2)
# (4) 定常：ZnO a, ZnSO4 b=40, 抜き取り V（Zn 1.0, H2SO4 1.0 mol/L）
a,V=symbols('a V'); b=40; zz=950
sol=solve([a+b-zz-1*V, b-2*V],[a,V],dict=True)[0]; print(sol)
acid=2*zz-2*sol[a]-2*1*sol[V]; print('H+ balance',acid)
for eq in ['2ZnS + 3O2 -> 2ZnO + 2SO2','ZnO + H2SO4 -> ZnSO4 + H2O','2ZnSO4 + 2H2O -> 2Zn + 2H2SO4 + O2','2H2O -> O2 + 4H+ + 4e-','Zn^2+ + 2e- -> Zn','2H+ + 2e- -> H2','H2SO4 + ZnO -> ZnSO4 + H2O']:
    print(eq,check_equation(eq))
