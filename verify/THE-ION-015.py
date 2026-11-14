import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, symbols, solve
from scipy.optimize import brentq
from math import log10
K1=R(5,1000); h=R(1,100)
Capp=h**2/K1; assert Capp==R(2,100) and h/Capp==R(1,2)
C=symbols('C'); Cex=solve(h**2/(C-h)-K1, C)[0]; print(Cex, float(h/Cex)); assert Cex==R(3,100)
k1,k2,k3,Kw=5e-3,6.3e-8,4.5e-13,1e-14; c=0.03
def f(lh):
    x=10**lh; D=x**3+k1*x**2+k1*k2*x+k1*k2*k3
    return x-Kw/x-c*(k1*x**2+2*k1*k2*x+3*k1*k2*k3)/D
print('exact pH', -brentq(f,-14,0)); assert round(-brentq(f,-14,0),2)==2.00
assert 3*Cex*10/(h*10)==9
assert not check_item('THE-ION-015')
