from sympy import symbols, solve, simplify, Rational as R
import numpy as np
from scipy.integrate import solve_ivp
E,EA,EB,A,B,ET,k1a,k2a,k3a,k1b,k2b,k3b=symbols('E EA EB A B ET k1a k2a k3a k1b k2b k3b',positive=True)
sol=solve([k1a*E*A-(k2a+k3a)*EA, k1b*E*B-(k2b+k3b)*EB, ET-E-EA-EB],[E,EA,EB],dict=True)[0]
KA=(k2a+k3a)/k1a; KB=(k2b+k3b)/k1b
vA=k2a*sol[EA]; vB=k2b*sol[EB]
assert simplify(vA - k2a*ET*A/(KA*(1+B/KB)+A))==0
assert simplify(vA/vB - (k2a/KA*A)/(k2b/KB*B))==0
def v_single(k2,K,S): return k2*S/(K+S)
a3=v_single(50,1e-4,1e-2); b3=v_single(400,2e-3,1e-2); print(a3,b3)
assert round(a3)==50 and round(b3,-1)==330
a4=50*1e-2/(1e-4*(1+1e-2/2e-3)+1e-2); b4=400*1e-2/(2e-3*(1+1e-2/1e-4)+1e-2); print(a4,b4,a4/b4)
assert round(a4)==47 and round(b4)==19 and abs(a4/b4-2.5)<1e-9
# ODE with elementary steps
p=dict(k1a=1.0e6,k2a=50,k3a=50,k1b=4.0e5,k2b=400,k3b=400)
assert abs((p['k2a']+p['k3a'])/p['k1a']-1e-4)<1e-15 and abs((p['k2b']+p['k3b'])/p['k1b']-2e-3)<1e-15
Et=1e-8
def f(t,y):
    e,ea,eb,a,b=y
    ra1=p['k1a']*e*a-p['k3a']*ea; rb1=p['k1b']*e*b-p['k3b']*eb
    return [-ra1-rb1+p['k2a']*ea+p['k2b']*eb, ra1-p['k2a']*ea, rb1-p['k2b']*eb, -ra1, -rb1]
s=solve_ivp(f,[0,4000],[Et,0,0,1e-2,1e-2],method='LSODA',rtol=1e-9,atol=1e-20,dense_output=True)
for t in [1,100,1000,3000]:
    e,ea,eb,a,b=s.sol(t)
    print(t, a, b, (p['k2a']*ea)/(p['k2b']*eb), (p['k2a']*ea/a)/(p['k2b']*eb/b), a/b)
    assert abs((p['k2a']*ea/a)/(p['k2b']*eb/b)-2.5)<0.01
e,ea,eb,a,b=s.sol(1); assert abs(p['k2a']*ea/Et-a4)<0.5 and abs(p['k2b']*eb/Et-b4)<0.5
ratios=[s.sol(t)[3]/s.sol(t)[4] for t in np.linspace(1,4000,50)]
assert all(x>y for x,y in zip(ratios,ratios[1:]))
print('THE-EQL-063 OK')
