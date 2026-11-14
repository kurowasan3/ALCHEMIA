import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, nsolve
A=[4.00,3.40,3.00,2.80]; B=[3.00,1.80,1.00,0.60]; C=[0,0.60,1.00,1.20]
for i in range(1,4):
    dA=A[0]-A[i]; dB=B[0]-B[i]; dC=C[i]; print(dA,dB,dC); assert abs(dB-2*dA)<1e-9 and abs(dC-dA)<1e-9
K = (1.20/2)/((2.80/2)*(0.60/2)**2); approx(K, 4.76, 0.002)
x = symbols('x')
def eqC(nA, nB, V):
    from scipy.optimize import brentq; f=lambda t:(t/V)/(((nA-t)/V)*((nB-2*t)/V)**2)-K; return brentq(f,1e-9,min(nA,nB/2)-1e-9)
c0 = eqC(4.0, 3.0, 2.0); print('C base', c0); approx(c0, 1.20, 0.01)
print('C at V=1.0', eqC(4.0, 3.0, 1.0)); assert eqC(4.0, 3.0, 1.0) > c0
print('C with A=5', eqC(5.0, 3.0, 2.0)); assert eqC(5.0, 3.0, 2.0) > c0
assert check_item('THE-EQL-019')
