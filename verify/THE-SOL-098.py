# 検証スクリプト（東工大 g12）  python3 -I この.py
import json, re
from sympy import symbols, Rational as R, solve, sqrt, nsimplify
ITEM = '/home/claude/chemx/tokodai/outh/items/THE-SOL-098/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
assert it['problem'].count('\\item') == len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])) == 4
S = R(5)            # cm^2
V0 = R(50)          # mL
RT = R(25,10)*10**6 # Pa L/mol
Pcm = 100           # Pa per cm
Kf = R(185,100); dT = R(740,10000); MNaCl = R(585,10)
# (1)
m_ion = dT/Kf; cNaCl = m_ion/2
nNaCl = cNaCl*R(100,1000)
mass = nNaCl*MNaCl; print('NaCl', float(mass)); assert abs(float(mass)-0.117) < 5e-4
# NaCl はつり合いの体積によらず 100 mL に均一：右管の濃度も同じ
# (2)
h = R(4)
VL = V0 + S*h/2; assert VL == 60
c = h*Pcm/RT; nP = c*VL/1000; M = R(120,1000)/nP
print('VL', VL, 'c', float(c), 'nP', float(nP), 'M', float(M)); assert M == 12500
wrong = R(120,1000)/(c*V0/1000); print('希釈無視', float(wrong)); assert wrong == 15000 and '1.50\\times10^{4}' in it['solution']
# (3)
x = symbols('x', positive=True)
sol = solve(Pcm*x - (2*nP)/((V0+S*x/2)/1000)*RT, x); print(sol, [float(s) for s in sol])
assert len(sol) == 1 and abs(float(sol[0]) - 7.09) < 5e-3
assert sol[0].equals(2*(sqrt(73)-5))
assert float(sol[0]) < 8
# (4) 水だけを通す膜
hmax = 2*V0/S; Pmax = hmax*Pcm; print('hmax', hmax, 'Pmax', Pmax); assert hmax == 20 and Pmax == 2000
Pi_min = (2*nNaCl + nP)/(R(100,1000))*RT; print('Π(最も薄い)', float(Pi_min)); assert Pi_min > Pmax
assert abs(float(Pi_min) - 1.0e5) < 1e3
# 途中のどの h（0〜20 cm）でも Π > 液柱の圧力
for k in range(0, 201):
    hh = R(k,10); V = V0 + S*hh/2
    assert (2*nNaCl+nP)/(V/1000)*RT > hh*Pcm
print('ALL OK')
