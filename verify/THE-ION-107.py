# THE-ION-107 検証：弱酸・強酸の pH の濃度依存（三次方程式で厳密解）
from sympy import *
h=symbols('h',positive=True)
Ka=Rational(1,10**5); Kw=Rational(1,10**14)
def pH_weak(C):
    sols=[s for s in Poly(h**3+Ka*h**2-(Kw+Ka*C)*h-Ka*Kw,h).nroots() if s.is_real and s>0]
    return -log(sols[0],10).evalf(), sols[0]
def pH_strong(C):
    H=(C+sqrt(C**2+4*Kw))/2; return -log(H,10).evalf()
for C in [Rational(1,10),Rational(1,10**3),Rational(1,10**5),Rational(1,10**8)]:
    p,H=pH_weak(C); print('C=',float(C),'weak pH',round(float(p),3),'alpha',float((H-Kw/H)/C),' strong pH',round(float(pH_strong(C)),3))
# 近似式
for C in [0.1,1e-3]:
    print('approx',0.5*(5-log(C,10)).evalf())
