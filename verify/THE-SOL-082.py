# THE-SOL-082 検証：冷却による析出量（水和物・溶解度が温度とともに減少する塩）
from sympy import symbols, solve, Rational as R
sol={'KNO3':(31.6,169.0),'Na2CO3':(21.5,44.0),'Li2CO3':(1.33,0.85)}  # (20C, 80C) g/100 g 水
add={'KNO3':40.0,'Na2CO3':30.0,'Li2CO3':1.2}
# KNO3
k=add['KNO3']-sol['KNO3'][0]; print('KNO3 dissolves at 80?', add['KNO3']<sol['KNO3'][1], 'ppt 20C', k); assert abs(k-8.4)<1e-9
# Na2CO3：析出は Na2CO3·10H2O（式量 106+180=286）
x=symbols('x')
r=solve((30-R(106,286)*x)/(100-R(180,286)*x)-R(215,1000),x)[0]
print('Na2CO3 dissolves at 80?', add['Na2CO3']<sol['Na2CO3'][1], 'decahydrate ppt', float(r)); assert round(float(r))==36
# 残った溶液の確認
rem_s=30-106/286*float(r); rem_w=100-180/286*float(r); print('check ratio', rem_s/rem_w*100)
# Li2CO3
u80=add['Li2CO3']-sol['Li2CO3'][1]; print('Li2CO3 undissolved at 80', u80, 'at 20 dissolves all?', add['Li2CO3']<=sol['Li2CO3'][0])
assert abs(u80-0.35)<1e-9 and add['Li2CO3']<=sol['Li2CO3'][0]
# Li2CO3 を 20 ℃で飽和させてから 80 ℃に加熱したときの析出量
print('Li2CO3 sat at 20 heated to 80: ppt', sol['Li2CO3'][0]-sol['Li2CO3'][1])
