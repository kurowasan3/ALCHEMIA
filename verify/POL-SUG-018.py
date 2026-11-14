import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it, b = common.basic("POL-SUG-018"); chk("basic", b)
# repeating unit = monomer - H2O
glcnac=CalcMolFormula(Chem.MolFromSmiles("OCC1OC(O)C(NC(C)=O)C(O)C1O")); glcn=CalcMolFormula(Chem.MolFromSmiles("OCC1OC(O)C(N)C(O)C1O"))
print(glcnac,glcn); chk("monomers C8H15NO6 / C6H13NO5 -> units C8H13NO5 / C6H11NO4", glcnac=="C8H15NO6" and glcn=="C6H13NO5")
chk("unit masses 203/161", 8*12+13+14+5*16==203 and 6*12+11+14+4*16==161)
x=symbols("x"); xs=solve(Q(1,10)*(203-42*x)-Q(1715,100),x)[0]; print(xs); chk("(1) 75%", xs==Q(3,4))
naoh=Q(1,10)*xs*40; chk("(2) 3.00 g", naoh==3)
print("RESULT:", "PASS" if ok_all else "FAIL")
