import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-PRT-014"); chk("basic", b)
K1,K2=Q(5,1000),Q(2,10**10)
r=(Q(1,10**4))**2/(K1*K2); print(r); chk("(2) 1.0e4 (cation>anion -> cathode)", r==10**4)
pI=-math.log10(math.sqrt(float(K1*K2))); print(pI); chk("(3) 6.0", round(pI,1)==6.0)
H=Q(1,10**9); a=K1/H; c=a*K2/H; print(a,c); chk("(4) 1:5.0e6:1.0e6", a==5*10**6 and c==10**6)
print("RESULT:", "PASS" if ok_all else "FAIL")
