# THE-ATM-033 verification: Rb-Sr isochron dating
import math
from sympy import symbols, Eq, solve, exp, log
lam=0.693/4.88e10
print('lambda',lam)
assert abs(lam-1.42e-11)<0.005e-11
R0=0.7050; t_true=1.20e9
f=math.exp(lam*t_true)-1
data={'A':0.500,'B':5.00,'C':12.0}
meas={k:round(R0+v*f,4) for k,v in data.items()}
print('measured (4 dp)',meas)
# solve from A and B
s=(meas['B']-meas['A'])/(data['B']-data['A'])
R0f=meas['A']-data['A']*s
lt=math.log(1+s); t=lt/1.42e-11
print('slope',s,'R0',R0f,'ln(1+s)',lt,'t',t,'approx t',s/1.42e-11)
assert round(s,4)==0.0172 and round(R0f,4)==0.7050 and round(t/1e9,2)==1.20
# check C prediction
predC=R0f+12.0*s; print('pred C',predC,meas['C'])
assert abs(predC-meas['C'])<2e-4
# mineral D (altered): rb/sr 8.00 measured 0.8100, predicted
predD=R0f+8.0*s; print('pred D',predD)
# young rock 1.0e6 y spread
sp=(12.0-0.5)*(math.exp(1.42e-11*1.0e6)-1); print('spread 1My',sp)
assert round(sp*1e4,1)==1.6
# exact symbolic check of isochron relation
N0,Np,l,tt=symbols('N0 Np l t',positive=True)
# parent now Np -> initial Np*exp(l t); daughter produced = Np*(exp(l t)-1)
print('OK')
