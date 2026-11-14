# POL-FIB-010 検証：界面重合（ナイロン610）
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
F = lambda s: re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles(s)))
assert F('ClC(=O)CCCCCCCCC(Cl)=O') == 'C10H16Cl2O2' and 10*12 + 16 + 71 + 32 == 239
assert F('NCCCCCCN') == 'C6H16N2' and 72 + 16 + 28 == 116
assert F('*C(=O)CCCCCCCCC(=O)NCCCCCCN*') == 'C16H30N2O2' and 16*12 + 30 + 28 + 32 == 282
na = R(239, 100) / 239; nb = R(870, 1000) / 116; print(na, nb)
lim = min(na, nb); m = 282 * lim; print(float(m)); assert lim == nb and round(float(m) + 1e-9, 2) == 2.12
assert 2 * lim < R(8, 10) / 40
assert not check_item('POL-FIB-010')
print('POL-FIB-010 OK')
