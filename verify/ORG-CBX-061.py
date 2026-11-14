# ORG-CBX-061 ヨードホルム反応でつながるアルコールと酸からなるエステル C9H18O2
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

for e in ['CH3CH2CH2CH(OH)CH3 + 4I2 + 6NaOH -> CHI3 + CH3CH2CH2COONa + 5NaI + 5H2O',
          '(CH3)2CHCH(OH)CH3 + 4I2 + 6NaOH -> CHI3 + (CH3)2CHCOONa + 5NaI + 5H2O',
          '(CH3)2CHCOOCH(CH3)CH(CH3)2 + H2O -> (CH3)2CHCOOH + (CH3)2CHCH(OH)CH3']:
    assert check_equation(e) is None, e

# R-COO-CH(CH3)-R で C 数 2r+3 = 9 → r = 3（C3H7 は 2 種）
esters = {'propyl': 'CCCC(=O)OC(C)CCC', 'isopropyl': 'CC(C)C(=O)OC(C)C(C)C'}
for k, s in esters.items():
    m = Chem.MolFromSmiles(s)
    f = CalcMolFormula(m)
    cc = FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    print(k, s, f, 'chiral centers', len(cc))
    assert f == 'C9H18O2' and len(cc) == 1
# r の他の値では C9 にならないことを確認
for r in range(1, 6):
    assert (2 * r + 3 == 9) == (r == 3)

# 脱水で生じるアルケン（立体異性体を区別）の数
def alkenes(alc):
    m = Chem.MolFromSmiles(alc)
    rw = Chem.RWMol(m)
    o = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c = m.GetAtomWithIdx(o).GetNeighbors()[0].GetIdx()
    prods = set()
    for nb in m.GetAtomWithIdx(c).GetNeighbors():
        if nb.GetSymbol() != 'C' or nb.GetTotalNumHs() == 0:
            continue
        r2 = Chem.RWMol(m)
        r2.GetBondBetweenAtoms(c, nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        r2.RemoveBond(c, o); r2.RemoveAtom(o)
        p = r2.GetMol(); Chem.SanitizeMol(p)
        for iso in EnumerateStereoisomers(p):
            prods.add(Chem.MolToSmiles(iso))
    return prods
p1 = alkenes('CCCC(O)C'); p2 = alkenes('CC(C)C(O)C')
print('2-pentanol ->', p1); print('3-methyl-2-butanol ->', p2)
assert len(p1) == 3 and len(p2) == 2
print('all OK')
