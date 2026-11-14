import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-016"); chk("basic", b)
n=symbols("n"); nC=Q(2288,100)/44; nH=2*Q(702,100)/18
sol=solve((6*n+10)*nH-(10*n+8)*nC,n); print(sol, N(nC), N(nH)); chk("(2) n=7", sol==[7])
Mw=12*52+78+16*35; print(Mw); chk("C52H78O35 = 1262, mass consistent", Mw==1262 and Q(1262,100)/Mw*52==nC and Q(1262,100)/Mw*78==nH)
chk("(1) composition independent of n", all(Q(6*12,162)==Q(6*k*12,162*k) for k in range(5,9)))
print("RESULT:", "PASS" if ok_all else "FAIL")
