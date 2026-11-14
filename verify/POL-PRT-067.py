# POL-PRT-067 検証：グルタミン酸のラクタム化（ピログルタミン酸）
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
glu = Chem.MolFromSmiles('OC(=O)CCC(N)C(=O)O')
P = Chem.MolFromSmiles('OC(=O)C1CCC(=O)N1')
lysL = Chem.MolFromSmiles('NC1CCCCNC1=O')
for n, m in [('Glu', glu), ('P', P), ('lysine lactam', lysL)]:
    print(n, CalcMolFormula(m), [len(r) for r in m.GetRingInfo().AtomRings()])
assert CalcMolFormula(P) == 'C5H7NO3'
assert [len(r) for r in P.GetRingInfo().AtomRings()] == [5]
assert [len(r) for r in lysL.GetRingInfo().AtomRings()] == [7]
# 一級アミン（ニンヒドリン陽性の目安）の有無
patt = Chem.MolFromSmarts('[NX3;H2][CX4]')
print('Glu primary amine', glu.HasSubstructMatch(patt), 'P primary amine', P.HasSubstructMatch(patt))
assert glu.HasSubstructMatch(patt) and not P.HasSubstructMatch(patt)
# 質量
M = lambda f: {'C': 12, 'H': 1.0, 'N': 14, 'O': 16}
Mglu = 5*12 + 9*1.0 + 14 + 4*16; MP = 5*12 + 7*1.0 + 14 + 3*16
print(Mglu, MP, 14.7/Mglu*MP)
assert abs(14.7/Mglu*MP - 12.9) < 0.005
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('C5H9NO4 -> C5H7NO3 + H2O') is None
print('OK')
