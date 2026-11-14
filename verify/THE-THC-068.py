# THE-THC-068 検証：直鎖アルカンの燃焼エンタルピーの CH2 あたりの増分とシクロアルカンのひずみ
# 生成エンタルピー（気体，kJ/mol，文献値）から燃焼エンタルピー（H2O は液体）を計算して問題の値と照合する
from sympy import Rational as R
CO2,H2O=R('-393.5'),R('-285.8')
dfH={'C2H6':(2,6,R('-84.0')),'C3H8':(3,8,R('-104.7')),'C4H10':(4,10,R('-125.6')),'C5H12':(5,12,R('-146.8')),
     'cC4H8':(4,8,R('27.7')),'cC5H10':(5,10,R('-77.0')),'cC6H12':(6,12,R('-123.4'))}
dc={k:c*CO2+h/2*H2O-f for k,(c,h,f) in dfH.items()}
for k,v in dc.items(): print(k, float(v))
given={'C2H6':-1560,'C3H8':-2219,'C4H10':-2877,'C5H12':-3536,'cC4H8':-2745,'cC5H10':-3320,'cC6H12':-3952}
for k in given: assert abs(float(dc[k])-given[k])<1.0, k
a=R(given['C5H12']-given['C2H6'],3)
print('CH2 あたり', float(a))
assert round(float(a))==-659
est={n:n*a for n in (4,5,6)}
strain={4:-(given['cC4H8']-est[4]),5:-(given['cC5H10']-est[5]),6:-(given['cC6H12']-est[6])}
for n in (4,5,6): print(n, float(est[n]), float(strain[n]), float(strain[n]/n))
assert round(float(est[4]))==-2635 and round(float(est[5]))==-3293 and round(float(est[6]))==-3952
assert round(float(strain[4]))==110 and round(float(strain[5]))==27 and abs(float(strain[6]))<1
assert round(float(strain[4]/4),1)==27.6 and round(float(strain[5]/5),1)==5.3
print('OK')
# 問題では増分を有効数字3桁 659 kJ/mol として使う
a2=659
s2={n:-(given[k])-n*a2 for n,k in ((4,'cC4H8'),(5,'cC5H10'),(6,'cC6H12'))}
print('659 を使ったひずみ', s2, {n:s2[n]/n for n in s2})
assert s2=={4:109,5:25,6:-2}
print('OK2')
