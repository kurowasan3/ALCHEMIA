# POL-RUB-008 検証：加硫
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
m = Chem.MolFromSmiles('*CC(C)=CC*')
assert re.sub(r'\*\d*', '', CalcMolFormula(m)) == 'C5H8' and sum(b.GetBondTypeAsDouble() == 2 for b in m.GetBonds()) == 1
db = R(68) / 68
for s, want in ((R(32, 10), 90), (R(192, 10), 40)):
    S = s / 32; cross = S / 2; used = 2 * cross
    rem = (db - used) / db * 100; print(s, rem); assert rem == want
assert not check_item('POL-RUB-008')
print('POL-RUB-008 OK')
