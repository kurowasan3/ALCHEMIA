# POL-PRT-021 検証：相補鎖・転写・水素結合数
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
s = 'TACGGATC'
comp = ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[b] for b in s)
rna = ''.join({'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}[b] for b in s)
print(comp, rna); assert comp == 'ATGCCTAG' and rna == 'AUGCCUAG'
def hb(frac, base):
    n = 200; x = round(n * frac); other = (n - 2 * x) // 2
    at, gc = (x, other) if base == 'A' else (other, x)
    return 2 * at + 3 * gc
print(hb(0.30, 'A'), hb(0.30, 'G')); assert (hb(0.30, 'A'), hb(0.30, 'G')) == (240, 260)
assert CalcMolFormula(Chem.MolFromSmiles('OCC(O)C(O)CC=O')) == 'C5H10O4'
assert CalcMolFormula(Chem.MolFromSmiles('OCC(O)C(O)C(O)C=O')) == 'C5H10O5'
assert not check_item('POL-PRT-021')
print('POL-PRT-021 OK')
