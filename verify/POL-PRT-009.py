import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it, b = common.basic("POL-PRT-009"); chk("basic", b)
F={k:CalcMolFormula(Chem.MolFromSmiles(s)) for k,s in {"Gly":"NCC(=O)O","GlyOEt":"NCC(=O)OCC","AcGly":"CC(=O)NCC(=O)O","propanoic":"CCC(=O)O"}.items()}
print(F); chk("formulas", F=={"Gly":"C2H5NO2","GlyOEt":"C4H9NO2","AcGly":"C4H7NO3","propanoic":"C3H6O2"})
chk("MW Gly 75 / propanoic 74", 2*12+5+14+32==75 and 3*12+6+32==74)
pat_amine=Chem.MolFromSmarts("[NX3;H2][CX4]"); pat_acid=Chem.MolFromSmarts("C(=O)[OH]")
fa={k:(Chem.MolFromSmiles(s).HasSubstructMatch(pat_amine),Chem.MolFromSmiles(s).HasSubstructMatch(pat_acid)) for k,s in {"Gly":"NCC(=O)O","GlyOEt":"NCC(=O)OCC","AcGly":"CC(=O)NCC(=O)O"}.items()}
print(fa); chk("(2) ester: amine only; AcGly: acid only", fa["GlyOEt"]==(True,False) and fa["AcGly"]==(False,True))
chk("(3) ninhydrin (free NH2): Gly, GlyOEt", [k for k,v in fa.items() if v[0]]==["Gly","GlyOEt"])
print("RESULT:", "PASS" if ok_all else "FAIL")
