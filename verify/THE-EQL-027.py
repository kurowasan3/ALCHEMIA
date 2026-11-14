import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item
from scipy.optimize import brentq
Kc = 5.0; N = 1.0  # NO2 換算の全量 (mol)
def eq(V):
    y = brentq(lambda y: (y/V)/(((N-2*y)/V)**2) - Kc, 0, N/2-1e-12); return (N-2*y), y  # n(NO2), n(N2O4)
nA, _ = eq(1.0); cA = nA/1.0
cB = nA/0.5
nC, _ = eq(0.5); cC = nC/0.5
nD, _ = eq(2.0); cD = nD/2.0
print(cA, cB, cC, cD, nA, nC)
assert cB > cC > cA and nA > nC and cD < cA and cD > cA/2
assert check_item('THE-EQL-027')
