# THE-SOL-091 検証：浸透圧から沸点上昇・凝固点降下
R = 8.31e3; T = 300
Pi = 2.49e5
c = Pi/(R*T); print('c', c)
assert round(c, 3) == 0.100
Kb, Kf = 0.52, 1.86
print('dTb', Kb*0.100, 'dTf', Kf*0.100)
assert round(Kb*0.100, 3) == 0.052 and round(Kf*0.100, 2) == 0.19
# NaCl：粒子濃度 0.100 → NaCl 0.0500 mol/L
print('NaCl', 0.100/2)
# タンパク質 1.0 g/L, M=5.0e4
cp = 1.0/5.0e4
print('protein Pi', cp*R*T, 'dTf', Kf*cp)
assert round(cp*R*T) == 50
print('OK')
