import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
it, b = common.basic("POL-SUG-015"); chk("basic", b)
def nch(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True))
chain="OCC(=O)C(O)C(O)C(O)CO"; fur="OCC1(O)OC(CO)C(O)C1O"; pyr="OCC1(O)OCC(O)C(O)C1O"
def ring(s):
    m=Chem.MolFromSmiles(s); return [len(r) for r in m.GetRingInfo().AtomRings()]
print(nch(chain),nch(fur),nch(pyr),ring(fur),ring(pyr))
chk("chain 3 chiral, furanose 4", nch(chain)==3 and nch(fur)==4)
chk("C5-O -> 5-ring, C6-O -> 6-ring", ring(fur)==[5] and ring(pyr)==[6])
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
chk("isomers C6H12O6", {CalcMolFormula(Chem.MolFromSmiles(s)) for s in (chain,fur,pyr)}=={"C6H12O6"})
print("RESULT:", "PASS" if ok_all else "FAIL")
