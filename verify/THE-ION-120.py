# THE-ION-120 verification: carrier equilibrium np = ni^2 with donors/acceptors
import sympy as sp, math
ni=1.0e10
def carriers(ND,NA,ni2):
    # neutrality n + NA = p + ND, p = ni2/n  ->  n^2 - (ND-NA) n - ni2 = 0 (exact, sympy)
    D=sp.Rational(int(ND))-sp.Rational(int(NA)); n=D/2+sp.sqrt(D**2/4+sp.Rational(int(ni2)))
    nn=float(n); return nn, float(ni2)/nn
n,p=carriers(1e15,0,ni**2); print('(1)',n,p); assert abs(n/1e15-1)<1e-6 and abs(p/1e5-1)<1e-4
n,p=carriers(1e10,0,ni**2); print('(2)',n,p); assert f'{n:.1e}'=='1.6e+10' and f'{p:.1e}'=='6.2e+09'
Eg=1.12; kB=8.62e-5
ND=1e15; ni2T=1.1*ND*0.10*ND
T=sp.symbols('T',positive=True)
Ts=sp.solve(sp.Eq(2*sp.log(ni)-Eg/kB*(1/T-sp.Rational(1,300)),sp.log(ni2T)),T)[0]; print('(3) T',Ts); assert round(float(Ts),-1) in (570.0,580.0)
print('   2 sig', f'{float(Ts):.1e}')
n,p=carriers(1e15,0.9e15,ni**2); print('(4a)',n,p,'ratio',(n+p)/(1e15+1e5)); assert abs((n+p)/(1e15)-0.10)<1e-3
n,p=carriers(1e15,1e15,ni**2); print('(4b)',n,p,'ratio',(n+p)/(1e15+1e5)); assert abs(n-1e10)<1 and abs((n+p)/1e15-2e-5)<1e-9
n,p=carriers(1e15,1.1e15,ni**2); print('(4c)',n,p); assert p>n and abs(p/1e14-1)<1e-6
print('OK')
