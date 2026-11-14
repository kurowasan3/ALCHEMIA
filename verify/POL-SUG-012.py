import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-012"); chk("basic", b)
from rdkit import Chem
M={"C":12,"H":1,"O":16}
def mw(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(M[a.GetSymbol()] for a in m.GetAtoms())
T=mw("COCC1OC(O)C(OC)C(OC)C1OC"); R=mw("COCC1OC(O)C(OC)C(OC)C1O"); D=mw("OCC1OC(O)C(OC)C(OC)C1O")
print(T,R,D); chk("MW 236/222/208", (T,R,D)==(236,222,208))
nT,nR,nD=Q(236,100)/236,Q(2664,100)/222,Q(208,100)/208
per=(nT+nR+nD)/nD; print(per); chk("(2) 14", per==14)
br=Q(227,100)*10**6/162/per; print(N(br)); chk("(3) 1.0e3", round(float(br),-2)==1000)
print("RESULT:", "PASS" if ok_all else "FAIL")
