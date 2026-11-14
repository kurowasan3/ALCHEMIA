# POL-PRT-075 verification: Boc/DCC dipeptide synthesis
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
eqs = ['C13H22N2 + H2O -> C13H24N2O',
       'C8H15NO4 + C3H7NO2 + C13H22N2 -> C11H20N2O5 + C13H24N2O',
       'C11H20N2O5 -> C6H12N2O3 + CO2 + C4H8',
       'C6H12N2O3 + H2O -> C5H10N2O3 + CH3OH']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
smi = {'Boc-Ala': ('CC(C)(C)OC(=O)NC(C)C(=O)O', 'C8H15NO4'),
       'H-Gly-OMe': ('NCC(=O)OC', 'C3H7NO2'),
       'DCC': ('C1CCC(CC1)N=C=NC1CCCCC1', 'C13H22N2'),
       'DCU': ('O=C(NC1CCCCC1)NC1CCCCC1', 'C13H24N2O'),
       'Boc-Ala-Gly-OMe': ('CC(C)(C)OC(=O)NC(C)C(=O)NCC(=O)OC', 'C11H20N2O5'),
       'H-Ala-Gly-OMe': ('NC(C)C(=O)NCC(=O)OC', 'C6H12N2O3'),
       'Ala-Gly': ('NC(C)C(=O)NCC(=O)O', 'C5H10N2O3'),
       'isobutene': ('C=C(C)C', 'C4H8')}
for k, (s, f) in smi.items():
    got = CalcMolFormula(Chem.MolFromSmiles(s)); print(k, got); assert got == f
M = lambda C, H, N, O: 12 * C + 1.0 * H + 14 * N + 16 * O
assert M(8, 15, 1, 4) == 189 and M(5, 10, 2, 3) == 146 and M(13, 24, 2, 1) == 224
n = 1.89 / 189
m1 = n * 0.80 * 0.90 * 0.95 * 146; m2 = n * 0.80 * 224
print(m1, m2); assert round(m1, 1) == 1.0 and round(m2, 1) == 1.8
print('OK')
