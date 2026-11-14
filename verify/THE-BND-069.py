# THE-BND-069 検証：岩塩型 KCl の格子面のイオンの面密度と面間隔（格子点の総当たり）
import itertools, math
a = 0.63
# 3x3x3 セルに拡張した格子点（分率座標）
Cl, K = [], []
R = range(-2, 4)
for i, j, k in itertools.product(R, R, R):
    for f in [(0,0,0),(.5,.5,0),(.5,0,.5),(0,.5,.5)]:
        Cl.append((i+f[0], j+f[1], k+f[2]))
    for f in [(.5,0,0),(0,.5,0),(0,0,.5),(.5,.5,.5)]:
        K.append((i+f[0], j+f[1], k+f[2]))

def density(points, normal, offset, region_test, area):
    n = 0
    for p in points:
        if abs(sum(pi*ni for pi, ni in zip(p, normal)) - offset) < 1e-9 and region_test(p):
            n += 1
    return n / area

# (100) 面 z=0，大きな領域 0<=x,y<L（端の重複を避けるため半開区間）
L = 2
in100 = lambda p: 0 <= p[0] < L and 0 <= p[1] < L
d100_Cl = density(Cl, (0,0,1), 0, in100, (L*a)**2)
d100_K = density(K, (0,0,1), 0, in100, (L*a)**2)
# (110) 面 x=y：面内座標 u=x（対角方向の長さ √2 a），z
in110 = lambda p: 0 <= p[0] < L and 0 <= p[2] < L
d110_Cl = density(Cl, (1,-1,0), 0, in110, (L*math.sqrt(2)*a)*(L*a))
d110_K = density(K, (1,-1,0), 0, in110, (L*math.sqrt(2)*a)*(L*a))
print('(100)', d100_Cl, d100_K, ' (110)', d110_Cl, d110_K)
assert round(d100_Cl, 1) == 5.0 and round(d100_K, 1) == 5.0
assert round(d110_Cl, 1) == 3.6 and round(d110_K, 1) == 3.6
# √2=1.41 を使った値
print('1.41/a^2 =', 1.41/a**2); assert round(1.41/a**2, 1) == 3.6
# (111) 面 x+y+z=1：Cl のみ
on111_Cl = [p for p in Cl if abs(sum(p)-1) < 1e-9]
on111_K = [p for p in K if abs(sum(p)-1) < 1e-9]
print('(111) Cl', len(on111_Cl) > 0, 'K', len(on111_K))
assert len(on111_K) == 0 and len(on111_Cl) > 0
# 面間隔：(100) 面は z = 0, 1/2, 1 ... ごとにイオンがある → a/2
zs = sorted(set(round(p[2], 3) for p in Cl + K if 0 <= p[2] <= 1))
print('z levels', zs); assert zs == [0, 0.5, 1]
d100 = a/2
# (110)：x-y = 0, ±1/2, ... → 間隔 (1/2)/√2 * a
vals = sorted(set(round(p[0]-p[1], 3) for p in Cl + K if -1 <= p[0]-p[1] <= 1))
print('x-y levels', vals)
d110 = 0.5 * a / math.sqrt(2)
print('d100 =', d100, ' d110 =', d110, ' (1.41 使用)', a*1.41/4)
assert round(d100, 2) == 0.32 and round(d110, 2) == 0.22 and round(a*1.41/4, 2) == 0.22
# 各 (110) 面も電気的に中性
assert d110_Cl == d110_K
print('OK')
