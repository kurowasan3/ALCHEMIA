"""ORG-HYD-017 検証：組成・オゾン分解"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from sympy import Rational, nsimplify
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, can
C, H = Rational(882, 10) / 12, Rational(118, 10) / 1
print('C:H =', float(C), float(H), float(H / C)); assert abs(float(H / C) - 8 / 5) < 0.01
print('C5H8 %C', float(Rational(60, 68) * 100), '%H', float(Rational(8, 68) * 100))
assert round(float(Rational(60, 68) * 100), 1) == 88.2 and round(float(Rational(8, 68) * 100), 1) == 11.8
rx = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.[C:2]=O')
rx_ring = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>([C:1]=O.[C:2]=O)')
def ozo_all(smi):
    """すべての C=C を切断した生成物の一覧"""
    m = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol(m)
    dbl = [b for b in rw.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
    for b in dbl:
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        rw.RemoveBond(i, j)
        for k in (i, j):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(k, o, Chem.BondType.DOUBLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm)
    return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mm, asMols=True))
for s in ['C1=CCCC1', 'C=CCC=C', 'C=CC=CC', 'C1=CCCCC1']:
    print(s, formula(s), '->', ozo_all(s), [formula(p) for p in ozo_all(s)])
assert ozo_all('C1=CCCC1') == [can('O=CCCCC=O')] and formula('O=CCCCC=O') == 'C5H8O2'
assert ozo_all('C=CCC=C').count('C=O') == 2 and can('O=CCC=O') in ozo_all('C=CCC=C')
assert len(set(ozo_all('C=CC=CC'))) == 3
assert ozo_all('C1=CCCCC1') == [can('O=CCCCCC=O')]
# 枝分かれのない鎖式 C5H8 ジエンのうち HCHO を 2 分子与えるもの
dienes = ['C=CC=CC', 'C=CCC=C', 'CC=C=CC', 'C=C=CCC']
ok = [d for d in dienes if '=C=' not in d and ozo_all(d).count('C=O') == 2]
print('HCHO x2:', ok); assert ok == ['C=CCC=C']
report_equations('ORG-HYD-017')
