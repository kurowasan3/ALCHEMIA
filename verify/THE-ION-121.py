# THE-ION-121 検証：チオ硫酸・アンモニアによるハロゲン化銀の溶解（配位子の消費を含む）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, nsimplify, sqrt, N as Nn
for e in ['AgBr + 2S2O3^2- -> [Ag(S2O3)2]^3- + Br-', 'AgCl + 2NH3 -> [Ag(NH3)2]+ + Cl-']:
    assert check_equation(e) is None
Ksp = {'AgCl': 1.8e-10, 'AgBr': 5.0e-13, 'AgI': 8.0e-17}
b_th, b_am = 3.0e13, 1.6e7
s = symbols('s', positive=True)
def solub(K, C0):
    # s^2 = K (C0 - 2s)^2,  0 < s < C0/2
    r = solve(s - sqrt(K)*(C0 - 2*s), s)
    return float(r[0])
res = {}
for X in Ksp:
    K = Ksp[X]*b_th
    ex = solub(K, 0.500)
    ap = (K**0.5)*0.500
    res[X] = ex
    print(X, 'K', K, 'sqrtK', K**0.5, 'exact s', ex, 'approx s', ap, 'free L', 0.5 - 2*ex, '[Ag+]', Ksp[X]/ex)
assert f"{Ksp['AgBr']*b_th:.2g}" == '15'
assert f"{res['AgBr']:.3g}" == '0.221'
assert f"{res['AgI']:.3g}" == '0.0223'
print('AgI approx', (Ksp['AgI']*b_th)**0.5*0.5, 'fraction ligand used', 2*res['AgI']/0.5)
# (4) フィルムの枚数
import math
nfilm = math.floor(res['AgBr']/1.5e-3)
print('films', res['AgBr']/1.5e-3, nfilm)
assert nfilm == 147
# (5) 6.0 mol/L アンモニア
for X in ['AgCl', 'AgBr']:
    K = Ksp[X]*b_am
    print('NH3', X, 'K', K, 'sqrtK', K**0.5, 's', solub(K, 6.0))
assert f"{solub(Ksp['AgCl']*b_am, 6.0):.2g}" == '0.29'
assert f"{solub(Ksp['AgBr']*b_am, 6.0):.2g}" == '0.017'
