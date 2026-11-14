import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item
from sympy import symbols, nsolve
# 2SO2 + O2 <=> 2SO3, Kp 一定（仮の値），SO2:O2=2:1 を 1 mol O2 基準で
x = symbols('x')
Kp = 5e-4  # 1/Pa（仮定）
def eq(P):
    s = float(nsolve((2*x)**2*(3-x)/((2-2*x)**2*(1-x)*P) - Kp, x, 0.5)); return s
for P in (1e5,):
    a = eq(P); b = eq(2*P); n1 = 3-a; n2 = 3-b
    vol_ratio = (n2/n1)/2; print('x', a, b, 'V ratio', vol_ratio, 'density ratio', 1/vol_ratio)
    assert 1/vol_ratio > 2
assert check_item('THE-EQL-022')
