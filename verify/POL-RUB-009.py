# POL-RUB-009 検証：ヒドロキシ酸とラクチド
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from g12lib import check_item
def chiral(s): return Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)
lac, C = 'CC(O)C(=O)O', 'OCCC(=O)O'
assert CalcMolFormula(Chem.MolFromSmiles(lac)) == CalcMolFormula(Chem.MolFromSmiles(C)) == 'C3H6O3'
assert len(chiral(lac)) == 1 and not chiral(C)
prim = Chem.MolFromSmarts('[CH2][OH]'); sec = Chem.MolFromSmarts('[CH1]([CH3])[OH]')
assert Chem.MolFromSmiles(C).HasSubstructMatch(prim) and Chem.MolFromSmiles(lac).HasSubstructMatch(sec)
def stereo(s):
    xs = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))}
    mir = lambda x: Chem.CanonSmiles(x.replace('@@', '!').replace('@', '@@').replace('!', '@'))
    return len(xs), sum(mir(x) == x for x in xs)
print('ラクチド', stereo('CC1OC(=O)C(C)OC1=O'), 'グリコリド', stereo('O=C1COC(=O)CO1'))
assert stereo('CC1OC(=O)C(C)OC1=O') == (3, 1) and stereo('O=C1COC(=O)CO1') == (1, 1)
assert CalcMolFormula(Chem.MolFromSmiles('CC1OC(=O)C(C)OC1=O')) == 'C6H8O4' and 72 + 8 + 64 == 144
n = R(180) / 90; assert n == 2 and 72 * n == 144 and 18 * n == 36
assert re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles('*OC(C)C(*)=O'))) == 'C3H4O2'
assert not check_item('POL-RUB-009')
print('POL-RUB-009 OK')
