import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-PRT-016"); chk("basic", b)
# 文献の pKa（代表値）で正味電荷を計算
pk={"Glu":dict(acid=[2.19,4.25],base=[9.67]),"Ala":dict(acid=[2.34],base=[9.69]),"Lys":dict(acid=[2.18],base=[8.95,10.53])}
def q(a,pH): return sum(1/(1+10**(pH-p)) for p in pk[a]["base"])-sum(1/(1+10**(p-pH)) for p in pk[a]["acid"])
def pI(a):
    lo,hi=0,14
    for _ in range(60):
        m=(lo+hi)/2
        if q(a,m)>0: lo=m
        else: hi=m
    return m
P={a:round(pI(a),1) for a in pk}; print(P, {a:round(q(a,6.0),2) for a in pk})
chk("pI order Glu<Ala<Lys (~3.2,6.0,9.7)", P["Glu"]<P["Ala"]<P["Lys"] and abs(P["Glu"]-3.2)<0.1 and abs(P["Ala"]-6.0)<0.1 and abs(P["Lys"]-9.7)<0.15)
chk("(1) pH6: Glu -, Ala ~0, Lys +", q("Glu",6)<-0.5 and abs(q("Ala",6))<0.05 and q("Lys",6)>0.5)
chk("pH2 all cationic, pH12 all anionic", all(q(a,2.0)>0 for a in pk) and all(q(a,12)<0 for a in pk))
order=sorted(pk,key=lambda a:P[a]); chk("(2) Glu,Ala,Lys / (3) reverse", order==["Glu","Ala","Lys"])
print("RESULT:", "PASS" if ok_all else "FAIL")
