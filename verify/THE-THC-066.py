# THE-THC-066 硝酸アンモニウム爆薬の酸素収支 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from fractions import Fraction as F
Hf={'NH4NO3':-366,'H2O':-242,'CO2':-394,'CO':-111,'N2':0,'O2':0,'C':0}
M={'NH4NO3':80,'C':12}
def dH(react,prod): return sum(n*Hf[s] for s,n in prod.items())-sum(n*Hf[s] for s,n in react.items())
cases={
 'AN only (per 2)': ({'NH4NO3':2},{'N2':2,'H2O':4,'O2':1}),
 'r=0.25 (4AN+C)': ({'NH4NO3':4,'C':1},{'N2':4,'H2O':8,'CO2':1,'O2':1}),
 'r=0.5 (2AN+C)':  ({'NH4NO3':2,'C':1},{'N2':2,'H2O':4,'CO2':1}),
 'r=1 (2AN+2C)':   ({'NH4NO3':2,'C':2},{'N2':2,'H2O':4,'CO':2}),
}
res={}
for k,(r,p) in cases.items():
    h=dH(r,p); m=sum(n*M[s] for s,n in r.items()); res[k]=(h,m,h/m*1000)
    print(k,h,m,round(h/m*1000))
assert res['AN only (per 2)'][0]==-236 and round(res['AN only (per 2)'][2])==-1475
assert res['r=0.5 (2AN+C)'][0]==-630 and round(res['r=0.5 (2AN+C)'][2])==-3663
assert round(res['r=0.25 (4AN+C)'][2])==-2608 and round(res['r=1 (2AN+2C)'][2])==-2489
# 気体の物質量 (1.00 kg の 2AN+C 混合物から, 水蒸気を含む)
n=7*1000/172; print('gas mol/kg',n); assert round(n,1)==40.7
for e in ['2NH4NO3 -> 2N2 + 4H2O + O2','2NH4NO3 + C -> 2N2 + 4H2O + CO2','4NH4NO3 + C -> 4N2 + 8H2O + CO2 + O2','NH4NO3 + C -> N2 + 2H2O + CO']:
    assert check_equation(e) is None, e
print('OK')
