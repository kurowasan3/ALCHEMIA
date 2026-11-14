from sympy import symbols, solve, simplify, Rational as R, sqrt, nsimplify
import numpy as np
S0,K1,CT,GC=symbols('S0 K1 CT GC',positive=True)
gc=solve(GC-K1*S0*(CT-GC),GC)[0]
assert simplify(gc-K1*S0*CT/(1+K1*S0))==0
# X
CTx=[0,2.0,4.0,6.0,8.0]; STx=[0.20,1.20,2.20,3.20,4.20]
K1v=5000; S0x=2.0e-4
for c,s in zip(CTx,STx):
    pred=(S0x+K1v*S0x/(1+K1v*S0x)*c*1e-3)*1e3; assert abs(pred-s)<1e-9
# Y
CTy=[0,1.2,2.8,4.8,7.2,10.0]; STy=[0.10,0.20,0.50,1.00,1.70,2.60]; S0y=0.10
r=[(s-S0y)/c for c,s in zip(CTy[1:],STy[1:])]; print([round(x,3) for x in r]); assert all(a<b for a,b in zip(r,r[1:]))
free=[c-2*(s-S0y) for c,s in zip(CTy,STy)]; print(free)
q=[(s-S0y)*1e-3/((f*1e-3)**2) for f,s in zip(free[1:],STy[1:])]; print(q); assert all(abs(x-100)<1e-6 for x in q)
K2=100/1.0e-4; assert abs(K2-1.0e6)<1e-3
# 1:1 fit residual for Y
A=np.vstack([np.array(CTy),np.ones(len(CTy))]).T; coef,res,_,_=np.linalg.lstsq(A,np.array(STy),rcond=None); print('1:1 linear fit residual',res)
assert res[0]>0.05
x=symbols('x',positive=True); xs=solve(200*x**2+x-R('0.020'),x)[0]; print(float(xs))
ST=100*xs**2+R('1.0e-4'); print(float(ST)); assert round(float(ST),4)==0.0062
print('THE-EQL-064 OK')
