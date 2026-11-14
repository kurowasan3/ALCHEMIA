import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item
# 数値例：P は Ea=80 kJ/mol, Q は Ea=40 kJ/mol，T0=400 K で交わる
R_ = 8.314; T0 = 400.0
def logk(Ea, T, lA): return -Ea/(2.303*R_*T) + lA
lAP = 5.0; lAQ = lAP - (80e3-40e3)/(2.303*R_*T0)
assert abs(logk(80e3, T0, lAP) - logk(40e3, T0, lAQ)) < 1e-9
rP = 10**(logk(80e3, 310, lAP) - logk(80e3, 300, lAP)); rQ = 10**(logk(40e3, 310, lAQ) - logk(40e3, 300, lAQ))
print('10K ratio P,Q', rP, rQ); assert rP > rQ
print('low T (350K) logk P,Q', logk(80e3, 350, lAP), logk(40e3, 350, lAQ)); assert logk(40e3, 350, lAQ) > logk(80e3, 350, lAP)
# 触媒：Ea 減少・A 同じ → どの T でも k 増加
assert all(logk(60e3, T, lAP) > logk(80e3, T, lAP) for T in (300, 400, 500))
assert check_item('THE-EQL-016')
