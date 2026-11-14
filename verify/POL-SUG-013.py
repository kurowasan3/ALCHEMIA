import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-013"); chk("basic", b)
x=symbols("x"); xs=solve(Q(1,10)*(162+42*x)-Q(2628,100),x)[0]; print(xs)
chk("(1) 2.40", xs==Q(12,5)); chk("(2) 80.0%", xs/3*100==80)
ac=Q(1,10)*xs*102; print(N(ac)); chk("(3) 24.5 g", round(float(ac),1)==24.5)
chk("unit masses 162 / +42", 6*12+10+5*16==162 and (2*12+3+16)-1==42)
print("RESULT:", "PASS" if ok_all else "FAIL")
