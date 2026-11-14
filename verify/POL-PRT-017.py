# POL-PRT-017 検証：トリペプチドの構成アミノ酸の決定
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g12')
from collections import Counter
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item

AA = {  # 略号: (SMILES, COOH数, フェノール性OH数, ベンゼン環)
 'Gly': ('NCC(=O)O', 1, 0, 0), 'Ala': ('NC(C)C(=O)O', 1, 0, 0), 'Ser': ('NC(CO)C(=O)O', 1, 0, 0),
 'Asp': ('NC(CC(=O)O)C(=O)O', 2, 0, 0), 'Glu': ('NC(CCC(=O)O)C(=O)O', 2, 0, 0),
 'Lys': ('NC(CCCCN)C(=O)O', 1, 0, 0), 'Phe': ('NC(Cc1ccccc1)C(=O)O', 1, 0, 1), 'Tyr': ('NC(Cc1ccc(O)cc1)C(=O)O', 1, 1, 1)}
# 問題文の分子式と RDKit の一致
given = {'Gly':'C2H5NO2','Ala':'C3H7NO2','Ser':'C3H7NO3','Asp':'C4H7NO4','Glu':'C5H9NO4','Lys':'C6H14N2O2','Phe':'C9H11NO2','Tyr':'C9H11NO3'}
for k,(smi,*_) in AA.items():
    assert CalcMolFormula(Chem.MolFromSmiles(smi)) == given[k], k

def comp(f):
    import re
    return Counter({e:int(n or 1) for e,n in re.findall(r'([A-Z][a-z]?)(\d*)', f)})

target = comp('C17H23N3O7') + comp('H4O2')
sols = []
for combo in itertools.combinations(AA, 3):
    tot = sum((comp(given[a]) for a in combo), Counter())
    cooh = sum(AA[a][1] for a in combo) - 2   # ペプチド結合に 2 個使われる
    phen = sum(AA[a][2] for a in combo)
    ring = sum(AA[a][3] for a in combo)
    if tot == target and cooh == 2 and cooh + phen == 2 and ring >= 1:
        sols.append(combo)
print('条件を満たす組合せ:', sols)
assert sols == [('Ser', 'Glu', 'Phe')]

# 配列（α結合のみ）の SMILES から分子式を確認
res = {'Ser': 'CO', 'Glu': 'CCC(=O)O', 'Phe': 'Cc1ccccc1'}
def pep(seq):
    s = 'N'
    for i, a in enumerate(seq):
        s += f'C({res[a]})C(=O)' + ('N' if i < len(seq) - 1 else 'O')
    return s
seqs = [p for p in itertools.permutations(['Ser', 'Glu', 'Phe']) if p[1] != 'Phe']
for p in seqs:
    f = CalcMolFormula(Chem.MolFromSmiles(pep(p)))
    print('-'.join(p), f)
    assert f == 'C17H23N3O7'
assert len(seqs) == 4
assert not check_item('POL-PRT-017')
print('POL-PRT-017 OK')
