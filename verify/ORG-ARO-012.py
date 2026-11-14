"""ORG-ARO-012 検証：組成式・凝固点降下・ジアゾ経路"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from sympy import Rational as Q
from g09lib import report_equations, formula
r = [Q(636, 10) / 12, Q(60, 10) / 1, Q(93, 10) / 14, Q(212, 10) / 16]
mn = min(r); ratio = [float(x / mn) for x in r]
print('ratio', [round(x, 2) for x in ratio]); assert [round(x) for x in ratio] == [8, 9, 1, 2]
M = Q(512, 100) * Q(755, 1000) / (Q(256, 1000) * Q(1, 10))
print('M', M, float(M)); assert round(float(M)) == 151
Mf = 8 * 12 + 9 + 14 + 32; print('C8H9NO2', Mf, [round(100 * x / Mf, 2) for x in (96, 9, 14, 32)])
assert Mf == 151
S = {'X': 'CCc1ccc([N+](=O)[O-])cc1', 'Y': 'CCc1ccc(N)cc1', 'W': 'CCc1ccc(N=Nc2ccc(O)cc2)cc1', 'phenol': 'CCc1ccc(O)cc1'}
for k, s in S.items(): print(k, formula(s))
assert formula(S['X']) == 'C8H9NO2' and formula(S['W']) == 'C14H14N2O' and formula(S['phenol']) == 'C8H10O'
report_equations('ORG-ARO-012')
