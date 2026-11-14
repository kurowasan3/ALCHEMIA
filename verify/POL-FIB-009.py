# POL-FIB-009 検証：ベンゼンからナイロン66
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
F = lambda s: re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles(s)))
exp = {'c1ccccc1': 'C6H6', 'C1CCCCC1': 'C6H12', 'OC1CCCCC1': 'C6H12O', 'O=C1CCCCC1': 'C6H10O', 'OC(=O)CCCCC(=O)O': 'C6H10O4',
       'N#CCCCCC#N': 'C6H8N2', 'NCCCCCCN': 'C6H16N2', '*NCCCCCCNC(=O)CCCCC(*)=O': 'C12H22N2O2'}
for s, f in exp.items():
    assert F(s) == f, (s, F(s))
assert 6 * 12 + 6 == 78 and 12 * 12 + 22 + 28 + 32 == 226
nb = R(156000, 78); m3 = nb / 2 * 226 / 1000; m4 = nb * 226 / 1000; print(m3, m4)
assert m3 == 226 and m4 == 452
assert not check_item('POL-FIB-009')
print('POL-FIB-009 OK')
