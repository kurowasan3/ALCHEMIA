import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-010"); chk("basic", b)
Cu2O=Q(635,10)*2+16; print("Cu2O",Cu2O)
a1=Q(342,100)/342*Cu2O; a2=2*Q(342,100)/342*Cu2O
nm=Q(572,1000)/Cu2O; pct=nm*342/Q(684,100)*100; a3=Q(684,100)/342*2*Cu2O
print(N(a1),N(a2),N(pct),N(a3))
chk("(1) 1.43", round(float(a1),2)==1.43); chk("(2) 2.86", round(float(a2),2)==2.86)
chk("(3) 20.0 %, 5.72 g", round(float(pct),1)==20.0 and round(float(a3),2)==5.72)
print("RESULT:", "PASS" if ok_all else "FAIL")
