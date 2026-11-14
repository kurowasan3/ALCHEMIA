# THE-ION-093 検証：HCrO4- と Cr2O7^2- の二量化平衡と希釈（K1=3.2e-7, K2=100）
from sympy import symbols, solve
K2=100.0; K1=3.2e-7
def frac(C,H=1e-3):
    y=symbols('y',positive=True)
    ys=[s for s in solve(2*K2*y**2+y*(1+K1/H)-C,y) if s>0][0]
    return float(ys), float(K2*ys**2), float(2*K2*ys**2/C)
def frac_simple(C):
    y=symbols('y',positive=True)
    ys=[s for s in solve(2*K2*y**2+y-C,y) if s>0][0]
    return float(ys), float(K2*ys**2), float(2*K2*ys**2/C)
print('(1)', K1/1e-3)
a=frac_simple(0.10); print('(2)',a); assert abs(a[0]-0.02)<1e-12 and abs(a[1]-0.04)<1e-12 and abs(a[2]-0.8)<1e-12
b=frac_simple(1e-3); print('(3)',b); assert round(b[2]*100)==15
c=frac_simple(0.20); print('C=0.20',c); assert round(c[2]*100)==85
for H in [1e-2,1e-3,1e-9]: print('full H=',H, frac(0.10,H))
print('pH9 ratio', K1/1e-9)
