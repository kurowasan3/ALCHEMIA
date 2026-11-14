# POL-RUB-010 検証：組成計算と，でたらめな配列・ブロック配列の鎖を実際に作ってオゾン分解生成物を数える
import random
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
M = {'H': 1.0, 'C': 12, 'O': 16}
def mw(f):
    import re
    return sum(M[e]*(int(n) if n else 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', f))
assert mw('C8H8') == 104 and mw('C4H6') == 54 and mw('CH2O') == 30 and mw('C4H6O2') == 86
assert CalcMolFormula(Chem.MolFromSmiles('O=CCCC=O')) == 'C4H6O2'
nBd = 1.68/22.4; s = (5.09 - nBd*54)/104; b = 0.450/30; a = nBd - b
print(nBd, s, a, b)
assert (round(s, 5), round(a, 5), round(b, 5)) == (0.01, 0.06, 0.015)
nP = 3.61/86; nQ = 4.64/86
print('succ P,Q', nP, nQ)
assert f'{nP:.3g}' == '0.042' and f'{nQ:.3g}' == '0.054'
LP = a/(a-0.0420); LQ = a/(a-0.0540); print('L', LP, LQ)
assert f'{LP:.3g}' == '3.33' and round(LQ) == 10
x = a/(a+b+s); print('x', x, 'random succ', a*x, 'L', 1/(1-x))
assert f'{a*x:.3g}' == '0.0424' and f'{1/(1-x):.3g}' == '3.4'

# 模擬：単位列を作り，オゾン分解で HCHO（1,2-単位の数）とスクシンアルデヒド（14-14 の隣接数）を数え，連鎖数＝単位数－隣接数を確認
def ozonolysis_counts(seq):
    hcho = seq.count('V')
    succ = sum(1 for i in range(len(seq)-1) if seq[i] == 'B' and seq[i+1] == 'B')
    runs = sum(1 for i in range(len(seq)) if seq[i] == 'B' and (i == 0 or seq[i-1] != 'B'))
    return hcho, succ, runs
random.seed(0)
N = 170000  # S:B:V = 2:12:3
seq = ['S']*20000 + ['B']*120000 + ['V']*30000
random.shuffle(seq)
h, sc, r = ozonolysis_counts(seq)
assert r == seq.count('B') - sc
print('random: succ/B', sc/120000, 'expected x', 12/17, 'mean run', 120000/r)
assert abs(sc/120000 - 12/17) < 0.01
# 1 つの B--B の具体例を SMILES でオゾン分解して OHC-CH2-CH2-CHO が出ることを確認
from rdkit.Chem import AllChem
rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.O=[C:2]')
mol = Chem.MolFromSmiles('CC/C=C/CC/C=C/CC')   # 2 つの 1,4-単位が隣接したモデル（両端はエチル）
m = Chem.RWMol(mol)
# 全 C=C を切る
bonds = [bd for bd in m.GetBonds() if bd.GetBondTypeAsDouble() == 2]
for bd in bonds:
    i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
    m.RemoveBond(i, j)
    for k in (i, j):
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(k, o, Chem.BondType.DOUBLE)
frags = sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(m.GetMol(), asMols=True, sanitizeFrags=True))
print(frags); assert 'O=CCCC=O' in frags
print('OK')
