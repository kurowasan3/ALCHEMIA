import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from sympy import nsolve, Symbol
it, b = common.basic("THE-ION-023"); chk("basic", b)
pK2=-math.log10(2.5e-10); print(pK2)
pH2=9.60+(0.30-0.48); pH3=9.60+0.48-0.85
chk("(2) 9.4 (given logs)", round(pH2,1)==9.4); chk("(3) 9.2 (given logs)", round(pH3,1)==9.2)
chk("exact logs agree", round(pK2+math.log10(4/6),1)==9.4 and round(pK2+math.log10(3/7),1)==9.2)
# 厳密計算（電荷収支，0.1 L, K1=4.5e-3, Kw=1e-14）
def exact(ctot,na,cl):
    K1,K2,Kw=4.5e-3,2.5e-10,1e-14
    def g(h):
        D=h*h/K1+h+K2
        return na+h+ctot*(h*h/K1)/D-cl-Kw/h-ctot*K2/D
    lo,hi=1e-14,1e-1   # g(lo)<0, g(hi)>0
    for _ in range(200):
        m=math.sqrt(lo*hi)
        if g(m)>0: hi=m
        else: lo=m
    return -math.log10(m)
e2=exact(0.10,0.040,0.0); e3=exact(0.10,0.040,0.010); print(e2,e3)
chk("exact equilibrium consistent", round(e2,1)==9.4 and round(e3,1)==9.2)
print("RESULT:", "PASS" if ok_all else "FAIL")
