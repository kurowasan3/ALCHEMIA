"""THE-THC-017 検証：共鳴エネルギーと結合エネルギー"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from sympy import Integer
from g09lib import report_equations
h1, h2, h3 = Integer(-120), Integer(-232), Integer(-208)
pred = 3 * h1; RE = h3 - pred
print('pred', pred, 'resonance', RE); assert pred == -360 and RE == 152
conj = h2 - 2 * h1; print('diene stabilization', conj); assert conj == 8
a = h3 - h2; b = h2 - h1; print('a', a, 'b', b); assert a == 24 and b == -112
assert a + b + h1 == h3   # ヘスの法則の整合
E = dict(CC=348, CdC=614, CH=413, BrBr=193, CBr=285, HBr=366)
sub = (E['CH'] + E['BrBr']) - (E['CBr'] + E['HBr'])
add_loc = (E['CdC'] + E['BrBr']) - (E['CC'] + 2 * E['CBr'])
add = add_loc + RE
print('substitution', sub, 'addition localized', add_loc, 'addition benzene', add)
assert sub == -45 and add_loc == -111 and add == 41
report_equations('THE-THC-017')
