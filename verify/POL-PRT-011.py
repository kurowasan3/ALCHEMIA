import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it, b = common.basic("POL-PRT-011"); chk("basic", b)
NQ=Q(4,1000)*14/Q(292,1000)*100; NR=Q(2,1000)*14/Q(330,1000)*100; print(N(NQ),N(NR))
chk("(1) 19.2 / 8.48", round(float(NQ),1)==19.2 and round(float(NR),2)==8.48)
M1=Q(292,1000)/Q(4,1000); chk("(2) 73 < 75", M1==73)
M2=Q(292,1000)/Q(2,1000); MR=Q(330,1000)/Q(2,1000); print(M2,MR)
lys=Chem.MolFromSmiles("NCCCCC(N)C(=O)O"); phe=Chem.MolFromSmiles("NC(Cc1ccccc1)C(=O)O")
am={"C":12,"H":1,"N":14,"O":16}
def mw(m): m=Chem.AddHs(m); return sum(am[a.GetSymbol()] for a in m.GetAtoms())
def cp(m): mm=Chem.AddHs(m); return 100*12*sum(a.GetSymbol()=="C" for a in mm.GetAtoms())/mw(m)
print(CalcMolFormula(lys),mw(lys),cp(lys),CalcMolFormula(phe),mw(phe),cp(phe))
chk("(3) Lys C6H14N2O2 146, C 49.3%", CalcMolFormula(lys)=="C6H14N2O2" and mw(lys)==M2 and round(cp(lys),1)==49.3)
chk("(4) Phe 165, C 65.5%", mw(phe)==MR and round(cp(phe),1)==65.5)
chk("each 1 chiral C", all(len(Chem.FindMolChiralCenters(m,includeUnassigned=True))==1 for m in (lys,phe)))
print("RESULT:", "PASS" if ok_all else "FAIL")
