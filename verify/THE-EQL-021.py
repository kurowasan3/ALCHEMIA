import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
R=8.3e3; T=1000; P=1.0e5
n = P*13.6/(R*T); x = n-0.100; print(n, x)
approx(0.100-x, 0.0361, 0.005); approx(2*x, 0.128, 0.005)
pCO = 2*x/n*P; pCO2 = (0.100-x)/n*P; Kp = pCO**2/pCO2; print(Kp); approx(Kp, 2.8e5, 0.02)
Kc = Kp/(R*T); print(Kc); approx(Kc, 3.3e-2, 0.02)
assert x < 0.100   # 黒鉛 0.100 mol は残る
# (5) 黒鉛 0.030 mol
nC=0.030; Q = (2*nC/(0.1+nC)*P)**2/((0.1-nC)/(0.1+nC)*P); print('Q', Q); assert Q < Kp
approx(Q, 4.0e4, 0.02)
assert check_item('THE-EQL-021')
