import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-014"); chk("basic", b)
ALL={1,2,3,4,6}
def methyl(used): return frozenset(ALL-used-{1})   # 1位は加水分解でOHにもどる
res={}
for k in [2,3,4,6]: res[f"1->{k}"]=sorted([tuple(sorted(methyl({1}))),tuple(sorted(methyl({k})))])
res["1<->1"]=sorted([tuple(sorted(methyl({1})))]*2)
for k,v in res.items(): print("  ",k,v)
X=sorted([(2,3,4,6),(2,3,4)]); Y=[(2,3,4,6),(2,3,4,6)]
chk("(1)(2) unique 1->6 for X", [k for k,v in res.items() if v==X]==["1->6"])
chk("(3) unique 1<->1 for Y", [k for k,v in res.items() if v==Y]==["1<->1"])
red={k:(k!="1<->1") for k in res}; chk("X reducing, Y not", red["1->6"] and not red["1<->1"])
print("RESULT:", "PASS" if ok_all else "FAIL")
