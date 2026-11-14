"""ORG-CBX-015 検証：収率計算"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from sympy import Rational as Q
from g09lib import report_equations
Ma = 7 * 12 + 6 + 32; Me = 8 * 12 + 8 + 32
assert (Ma, Me) == (122, 136)
n0 = Q(122, 10) / Ma; theo = n0 * Me
y1 = Q(816, 100) / theo * 100
nr = Q(244, 100) / Ma; theo2 = (n0 - nr) * Me
y2 = Q(816, 100) / theo2 * 100
print('n0', n0, 'theo', float(theo), 'yield', float(y1), '| recovered', nr, 'theo2', float(theo2), 'yield2', float(y2))
assert y1 == 60 and y2 == 75
report_equations('ORG-CBX-015')
