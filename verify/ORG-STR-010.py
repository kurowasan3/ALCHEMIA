import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem import AllChem
def stereo(smi): return sorted({Chem.MolToSmiles(m) for m in EnumerateStereoisomers(Chem.MolFromSmiles(smi))})
def mirror(smi): return Chem.MolToSmiles(Chem.MolFromSmiles(smi.replace("@@","!").replace("@","@@").replace("!","@")))
it, b = common.basic("ORG-STR-010"); chk("basic", b)
ald=stereo("OCC(O)C(O)C=O"); print(ald)
nchi=len(Chem.FindMolChiralCenters(Chem.MolFromSmiles("OCC(O)C(O)C=O"),includeUnassigned=True))
chk("(1) 2 chiral C, 4 stereoisomers", nchi==2 and len(ald)==4)
alc=stereo("OCC(O)C(O)CO"); print(alc); chk("(2) 3 stereoisomers", len(alc)==3)
# reduce each aldotetrose: C=O -> CH2OH with stereo kept
rxn=AllChem.ReactionFromSmarts("[CH1:1]=[O:2]>>[CH2:1][O:2]")
prods={}
for s in ald:
    p=rxn.RunReactants((Chem.MolFromSmiles(s),))[0][0]; Chem.SanitizeMol(p); ps=Chem.MolToSmiles(p)
    prods[s]=ps; print("  ",s,"->",ps,"meso" if mirror(ps)==ps else "chiral")
inactive=[s for s,p in prods.items() if mirror(p)==p]
chk("(3) 2 give inactive, 3 distinct products", len(inactive)==2 and len(set(prods.values()))==3)
g=Chem.MolFromSmiles("OCC1OC(O)C(O)C(O)C1O"); ch=Chem.MolFromSmiles("OCC(O)C(O)C(O)C(O)C=O")
n1=len(Chem.FindMolChiralCenters(g,includeUnassigned=True)); n0=len(Chem.FindMolChiralCenters(ch,includeUnassigned=True))
print(n0,n1); chk("(4) glucose chain 4, ring 5", n0==4 and n1==5)
print("RESULT:", "PASS" if ok_all else "FAIL")
