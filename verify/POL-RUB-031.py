# POL-RUB-031 検証：破過曲線の面積と交換容量，電荷のつり合い
from sympy import Rational as R, integrate, symbols, Piecewise
V=symbols('V')
c0=R(1,100)
cout=Piecewise((0,V<150),(c0*(V-150)/100,V<250),(c0,True))
held=integrate((c0-cout)/1000,(V,0,300))   # mol（V は mL）
print('保持された Mg2+ mol', held)
assert held==R(2,1000)
print('交換容量（Na+換算）', 2*held); assert 2*held==R(4,1000)
# 100 mL の時点：Mg2+ はすべて保持され，出口は Na+ 2*c0
assert 2*c0==R(2,100)
# 150〜250 mL の Na+ 濃度：2*(c0-cout) → 線形に 0.020 → 0
for v in (150,200,250): print(v, 2*(c0-cout.subs(V,v)))
# 流出した Na+ の総量
naout=integrate(2*(c0-cout)/1000,(V,0,300)); assert naout==R(4,1000)
# 濃度2倍：保持量は同じ，破過までの体積は半分
c1=2*c0; scale=c0/c1
print('2倍濃度で出始める体積', 150*scale, '飽和', 250*scale)
cout2=Piecewise((0,V<75),(c1*(V-75)/50,V<125),(c1,True))
assert integrate((c1-cout2)/1000,(V,0,300))==R(2,1000)
print('OK')
