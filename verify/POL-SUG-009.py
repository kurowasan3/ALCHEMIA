import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it, b = common.basic("POL-SUG-009"); chk("basic", b)
S={"マルトース":"OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O",
   "イソマルトース":"OCC1OC(OCC2OC(O)C(O)C(O)C2O)C(O)C(O)C1O",
   "トレハロース":"OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O",
   "スクロース":"OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
   "パラチノース":"OCC1OC(OCC2OC(O)(CO)C(O)C2O)C(O)C(O)C1O"}
pat=Chem.MolFromSmarts("[OX2H1;!R]-[C;R]-[O;R]")
for k,s in S.items():
    m=Chem.MolFromSmiles(s); n=len(m.GetSubstructMatches(pat)); print(" ",k,CalcMolFormula(m),"hemiacetal OH =",n)
    S[k]=(CalcMolFormula(m),n)
chk("all C12H22O11", all(v[0]=="C12H22O11" for v in S.values()))
chk("(1) reducing = maltose, isomaltose", [k for k in ["マルトース","イソマルトース","トレハロース"] if S[k][1]>0]==["マルトース","イソマルトース"])
chk("(2) MeOH: maltose 1, trehalose 0", S["マルトース"][1]==1 and S["トレハロース"][1]==0)
chk("(3) palatinose reducing, sucrose not", S["パラチノース"][1]==1 and S["スクロース"][1]==0)
print("RESULT:", "PASS" if ok_all else "FAIL")
