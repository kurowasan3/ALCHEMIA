# THE-GAS-043 検証：密閉容器（定積）中の水の加熱（平坦部のない加熱曲線）
from sympy import nsolve, symbols, Rational as Q
R = 8.31e3; V = 10.0
ps = {40:7.38e3, 50:12.35e3, 60:19.95e3, 70:31.20e3, 80:47.41e3, 90:70.18e3, 100:101.3e3}
TK = lambda t: t + 273
cap = lambda t: ps[t]*V/(R*TK(t))          # 気体として存在できる最大の物質量
for t in ps: print(t, 'n_max =', round(cap(t),4))
n = 0.100
# (1) 液体が消える温度
lo = max(t for t in ps if cap(t) < n); hi = min(t for t in ps if cap(t) >= n)
print('(1) between', lo, hi)
assert (lo, hi) == (60, 70)
f = lambda t: ps[t] - n*R*TK(t)/V
tx = lo + 10*(-f(lo))/(f(hi)-f(lo)); print('   nRT/V 60,70:', round(n*R*TK(60)/V/1e3,3), round(n*R*TK(70)/V/1e3,3), ' 差', round(f(60)/1e3,3), round(f(70)/1e3,3)); print('   直線補間 t =', round(tx,2))
assert 67.0 < tx < 67.6
# エネルギーのモデル（40 ℃ 液体基準）
Cl, Cg, q40 = 75.0, 25.0, 4.10e4
def E(t, ntot):
    ng = min(cap(t), ntot); nl = ntot - ng
    return nl*Cl*(t-40) + ng*(q40 + Cg*(t-40))
Q4060 = E(60,n)-E(40,n); print('(2) Q(40->60) =', round(Q4060,1), 'J ; 蒸発なし', n*Cl*20)
assert abs(Q4060-1.9e3) < 60
print('    Q(40->50) =', round(E(50,n)-E(40,n),1), ' Q(50->60) =', round(E(60,n)-E(50,n),1))
assert E(60,n)-E(50,n) > E(50,n)-E(40,n)          # 傾きは次第に小さくなる
Q7090 = E(90,n)-E(70,n); print('(3) Q(70->90) =', Q7090)
assert abs(Q7090 - 50) < 1e-9
print('    1 K あたりの熱の比 =', round((Q4060/20)/(Q7090/20),1))
# (4) 80 ℃ で折れる水の量と 100 ℃ で液体が残る条件
n4 = cap(80); print('(4) n =', round(n4,4), 'mol =', round(n4*18,3), 'g')
n100 = cap(100); print('    100 ℃ で液体が残る: n >', round(n100,4), 'mol =', round(n100*18,2), 'g')
assert round(n4*18,1) == 2.9 and round(n100*18,1) == 5.9
print('OK')
