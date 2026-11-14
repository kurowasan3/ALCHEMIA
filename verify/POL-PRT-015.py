import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from chemlib import check_equation
it, b = common.basic("POL-PRT-015"); chk("basic", b)
am={"C":12,"H":1.0,"N":14,"O":16}
E=Chem.MolFromSmiles("NCC(=O)NC(Cc1ccccc1)C(=O)OC"); f=CalcMolFormula(E); print(f)
cnt={"C":12,"H":16,"N":2,"O":3}; M=sum(am[k]*v for k,v in cnt.items())
pc={k:round(am[k]*v/M*100,1) for k,v in cnt.items()}; print(M,pc)
chk("E C12H16N2O3, M 236, % matches", f=="C12H16N2O3" and M==236 and pc=={"C":61.0,"H":6.8,"N":11.9,"O":20.3})
r=[61.0/12,6.8/1.0,11.9/14,20.3/16]; mn=min(r); rat=[x/mn for x in r]; print([round(x*2,2) for x in rat])
chk("empirical ratio 12:16:2:3", [round(x*2) for x in rat]==[12,16,2,3])
eq=check_equation("C12H16N2O3 + 2H2O -> C2H5NO2 + C9H11NO2 + CH4O"); print(" eq(B=C9H11NO2):",eq); chk("hydrolysis balanced", eq is None)
D=CalcMolFormula(Chem.MolFromSmiles("NC(Cc1ccccc1)C(=O)OC")); chk("D = Phe-OMe C10H13NO2", D=="C10H13NO2")
chk("(4) Gly-OMe C3H7NO2", CalcMolFormula(Chem.MolFromSmiles("NCC(=O)OC"))=="C3H7NO2")
chk("isomer E' same formula", CalcMolFormula(Chem.MolFromSmiles("NC(Cc1ccccc1)C(=O)NCC(=O)OC"))=="C12H16N2O3")
print("RESULT:", "PASS" if ok_all else "FAIL")
