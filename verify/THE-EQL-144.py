# THE-EQL-144 二本鎖形成の平衡と融解温度の濃度依存（二分子会合と分子内ヘアピン）
from sympy import symbols, solve, Rational as R, log, nsimplify
import math
Rg=8.31; dH=-3.6e5; dS=-1.00e3; ln2=0.693; ln10=2.30
c0=symbols('c0',positive=True)
# Tm での条件：F+G<=>D, 各鎖の全濃度 c0, 半分が二本鎖
D=c0/2; F=c0-D; K=D/(F*F); print('K at Tm =',K)
def Tm(c,exact=False):
    lnK=math.log(2/c) if exact else ln2+(-math.log10(c))*ln10
    return dH/(dS-Rg*lnK)
for c in [1e-6,1e-4]:
    print(c,'Tm(approx ln)=%.2f exact=%.2f'%(Tm(c),Tm(c,True)))
print('hairpin Tm',dH/dS)
# (4) 逆算
T1,T2=320.0,330.0
dH4=-Rg*(2*ln10)/(1/T1-1/T2); print('dH from Tm pair %.4g J/mol'%dH4, 'exact ln',-Rg*math.log(100)/(1/T1-1/T2))
dS4=dH4/T1+Rg*(ln2+6*ln10); print('dS %.4g'%dS4)
# 検算: 1/Tm 線形で exact 解と一致するか（数値で平衡を解く）
x=symbols('x',positive=True)
for c,T in [(1e-6,T1),(1e-4,T2)]:
    Kv=math.exp(-dH4/(Rg*T)+dS4/Rg)
    s=[float(r) for r in solve(Kv*(c-x)**2-x,x) if 0<float(r)<c][0]
    print(c,T,'fraction paired',s/c)
