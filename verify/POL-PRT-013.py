import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it, b = common.basic("POL-PRT-013"); chk("basic", b)
AA={"Gly":dict(chiral=0,ring=0,S=0),"Ala":dict(chiral=1,ring=0,S=0),"Ser":dict(chiral=1,ring=0,S=0),
    "Cys":dict(chiral=1,ring=0,S=1),"Phe":dict(chiral=1,ring=1,S=0),"Glu":dict(chiral=1,ring=0,S=0)}
sols=[]
for t in itertools.permutations(AA,3):
    if not any(AA[a]["chiral"]==0 for a in t): continue
    d1,d2=t[:2],t[1:]
    def Ylike(d): return any(AA[a]["ring"] for a in d) and not any(AA[a]["S"] for a in d)
    def Zlike(d): return any(AA[a]["S"] for a in d) and not any(AA[a]["ring"] for a in d)
    if not ((Ylike(d1) and Zlike(d2)) or (Ylike(d2) and Zlike(d1))): continue
    if t[0]!="Phe": continue
    sols.append(t)
print(sols); chk("unique sequence Phe-Gly-Cys", sols==[("Phe","Gly","Cys")])
# without N-terminal info: 2 candidates
X=Chem.MolFromSmiles("NC(Cc1ccccc1)C(=O)NCC(=O)NC(CS)C(=O)O"); print(CalcMolFormula(X))
chk("formula C14H19N3O4S", CalcMolFormula(X)=="C14H19N3O4S")
pb=len(X.GetSubstructMatches(Chem.MolFromSmarts("C(=O)N[#6]"))); chk("2 peptide bonds (biuret +), dipeptide 1 (−)", pb==2)
chk("phenyllactic acid = Phe with NH2->OH", CalcMolFormula(Chem.MolFromSmiles("OC(Cc1ccccc1)C(=O)O"))=="C9H10O3")
print("RESULT:", "PASS" if ok_all else "FAIL")
