# THE-THC-057 断熱容器での吸熱的な溶解が飽和で止まる温度から溶解エンタルピーを求める
from sympy import Rational as R, N
s40, slope = R(64), R(2)          # 40.0 ℃ で 64.0 g/水100 g, 1 ℃ 下がるごとに 2.00 g 減少
C = R(640)                        # J/K
M = R(85)
T = R(25)
s25 = s40 - slope*(40 - T)
print('dissolved', s25); assert s25 == 34
assert 60 > s25                   # 固体が残る（加えた 60.0 g より少ない）
q = C*(40 - T)
dH = q/(s25/M)/1000
print('dH =', N(dH)); assert dH == 24
# (4) 溶解度が温度とともに増加 ⇔ 吸熱（ルシャトリエ）: 符号の整合
assert slope > 0 and dH > 0
print('all OK')
