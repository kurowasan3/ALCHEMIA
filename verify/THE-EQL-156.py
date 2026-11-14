# THE-EQL-156 検証：エチレンの水和の ΔG°・平衡定数・反応の向き・平衡転化率と凝縮の判定
import sys, math
sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, nsolve
print(check_equation('C2H4 + H2O -> C2H5OH') or 'OK')
R = 8.31; dH = -45.6e3; dS = -126.0
# 熱力学データからの確認（標準生成エンタルピー・標準エントロピーの文献値）
print('文献値 ΔH', -234.8-(52.4-241.8), 'ΔS', 281.6-(219.3+188.8))
T1 = dH/dS; print('(1) K=1 の温度', T1)
def K(T):
    dG = dH - T*dS
    return dG, 10**(-dG/(2.303*R*T))
for T in (573, 473):
    dG, k = K(T); print(T, 'ΔG°', dG, 'log K', -dG/(2.303*R*T), 'K', k)
# (3) 573 K で P(C2H4)=3.0 MPa, P(H2O)=3.0 MPa, P(EtOH)=0.050 MPa, P°=0.10 MPa
dG573, K573 = K(573)
Q = (0.050/0.10)/((3.0/0.10)*(3.0/0.10)); print('Q', Q, 'Q/K', Q/K573)
dG = dG573 + 2.303*R*573*math.log10(Q); print('ΔG', dG)
# (4) 全圧 6.0 MPa 一定，C2H4:H2O = 1:1 から平衡
x = symbols('x', positive=True)
def conv(Kv, P=60.0):
    s = [v for v in solve(x*(2-x)/(1-x)**2 - Kv*P, x) if 0 < v < 1]; return float(s[0])
x573 = conv(K573); print('x573', x573, '分圧 MPa', [6*(1-x573)/(2-x573), 6*(1-x573)/(2-x573), 6*x573/(2-x573)])
# 有効数字を丸めた K で
print('x573 (K=3.75e-3)', conv(3.75e-3))
dG473, K473 = K(473)
x473 = conv(K473); pw = 6*(1-x473)/(2-x473); pe = 6*x473/(2-x473)
print('x473', x473, 'P_H2O', pw, 'P_EtOH', pe)
# 水の飽和蒸気圧 473 K: 1.55 MPa, 573 K: 8.59 MPa；エタノール 473 K: 3.0 MPa
print('凝縮判定 473K 水', pw > 1.55, '初期の水の分圧 3.0 MPa >1.55', 'エタノール', pe > 3.0)
print('573K 水 3.0 < 8.59 で凝縮しない')
# (5) 473 K で全圧を下げて水が凝縮しない上限：平衡時の P_H2O = 6*(1-x)/(2-x)*(P/6) <= 1.55, 初期 P/2 <= 1.55 → P <= 3.1 MPa
Pmax = 3.1
xP = conv(K473, P=Pmax/0.10); print('473K, 3.1 MPa の x', xP, 'P_H2O 平衡', Pmax*(1-xP)/(2-xP))
