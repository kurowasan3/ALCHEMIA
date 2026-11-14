# THE-EQL-162 verification: Curtin-Hammett / Halpern-type selectivity, steady-state formula and full ODE check
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
k1,k2,km1,km2,x,y=sp.symbols('k1 k2 km1 km2 x y',positive=True)
C,S=sp.symbols('C S',positive=True)
X=k1*C*S/(km1+x); Y=k2*C*S/(km2+y)
ratio=sp.simplify(x*X/(y*Y)); print('ratio',ratio)
r=sp.simplify(ratio.subs({k1:1,k2:1,km1:20,km2:1,y:x/1000})); print('P/Q =',sp.factor(r))
assert sp.simplify(r-(x+1000)/(x+20))==0
print('equilibrium Y/X', (1/1)/(1/20))  # K_Y/K_X = (k2/km2)/(k1/km1) = 20
low=sp.limit(r,x,0); high=sp.limit(r,x,sp.oo); print('low',low,'high',high); assert low==50 and high==1
ee=lambda q:(q-1)/(q+1)*100
print('ee low',float(ee(low)))
e2=float(ee(r.subs(x,2))); e20=float(ee(r.subs(x,20))); print('ee x=2',e2,'x=20',e20)
assert round(e2)==96 and round(e20)==92
xs=sp.solve(sp.Eq(r,19),x)[0]; print('x for ee90',xs,float(xs),'pressure Pa',float(xs)/2*1.0e5)
assert round(float(xs)/2*1.0e5,-5)==1.7e6
# full ODE simulation (catalyst 1e-3, S in excess, kon chosen so that binding is not saturating)
def sim(xv):
    kon=1.0; Ct=1e-3; S0=1.0
    def f(t,u):
        c,s,X_,Y_,P,Q=u
        v1=kon*c*s; v2=kon*c*s; r1=20*X_; r2=1.0*Y_; a=xv*X_; b=xv/1000*Y_
        return [-v1-v2+r1+r2+a+b, -v1-v2+r1+r2, v1-r1-a, v2-r2-b, a, b]
    sol=solve_ivp(f,(0,4000),[Ct,S0,0,0,0,0],method='LSODA',rtol=1e-10,atol=1e-14)
    P,Q=sol.y[4,-1],sol.y[5,-1]; return P/Q, P+Q
for xv in (0.01,2,20,34.44,1e4):
    q,tot=sim(xv); print('x',xv,'sim P/Q',q,'formula',float(r.subs(x,xv)),'conv',tot)
    assert abs(q/float(r.subs(x,xv))-1)<0.03
print('OK')
