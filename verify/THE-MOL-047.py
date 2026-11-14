# THE-MOL-047 多段の工程を通した炭素と酸素の原子の収支（銅の酸化物の還元）
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R
C=12;O=16;Si=28;Cu=64
coke=600; left=168; Si_g=R(840,10); O2=R(448)/R(224,10); ore=2420
nC=R(coke,C); nSi=Si_g/Si; nCc=2*nSi
nCa=nC-nCc-R(left,C); print('C total',nC,'C in C',nCc,'C in A',nCa)
x,y=symbols('x y')  # CO2, CO
s=solve([x+y-nCa,2*x+y-2*O2],[x,y]); print('A gas',s)
nO_removed=s[y]; mSiO2=nSi*(Si+2*O); mox=ore-mSiO2; mCu=mox-nO_removed*O
print('SiO2 g',mSiO2,'oxide g',mox,'Cu g',mCu,'Cu mol',mCu/Cu)
a,b=symbols('a b'); print(solve([a+2*b-mCu/Cu, a+b-nO_removed],[a,b]))
for e in ['2C + O2 -> 2CO','C + O2 -> CO2','CuO + CO -> Cu + CO2','Cu2O + CO -> 2Cu + CO2','SiO2 + 2C -> Si + 2CO']:
    print(e,check_equation(e))
