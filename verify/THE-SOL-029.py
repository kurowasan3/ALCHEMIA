import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("THE-SOL-029"); chk("basic", b)
M=Q(1)*Q(831,100)*10**3*300/(Q(247,100)*10**3*Q(2,10)); print(N(M))
chk("(1) 5.0e3", round(float(M),-2)==5000)
n=(M-18)/162-1; print(N(n)); chk("(2) n=30", round(float(n))==30)
chk("reverse: n=30 -> M=5040 -> Pi", abs(float(Q(1)/5040/Q(2,10)*Q(831,10)*100*300)-2470)<10)
chk("(3) 31x", 30+1==31)
print("RESULT:", "PASS" if ok_all else "FAIL")
