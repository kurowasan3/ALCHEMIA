# ORG-STR-058 検証：C5H10O のカルボニル化合物の非等価な炭素の数と反応による判別
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
def nC(smi):
    m = Chem.MolFromSmiles(smi); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol() == 'C'})
carb = {'ペンタナール': 'CCCCC=O', '2-メチルブタナール': 'CCC(C)C=O', '3-メチルブタナール': 'CC(C)CC=O',
        '2,2-ジメチルプロパナール': 'CC(C)(C)C=O', '2-ペンタノン': 'CCCC(C)=O', '3-ペンタノン': 'CCC(=O)CC',
        '3-メチル-2-ブタノン': 'CC(C)C(C)=O'}
# 網羅の確認：C5 骨格の C=O を全列挙（C4 アルデヒド骨格と C=O の位置）
from itertools import product
allc = set()
for sk in ['CCCCC', 'CC(C)CC', 'CC(C)(C)C']:
    m = Chem.MolFromSmiles(sk)
    for a in m.GetAtoms():
        if a.GetTotalNumHs() >= 2:
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(8)); rw.AddBond(a.GetIdx(), i, Chem.BondType.DOUBLE)
            try:
                mm = rw.GetMol(); Chem.SanitizeMol(mm); allc.add(Chem.MolToSmiles(mm))
            except Exception: pass
canon = {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in carb.values()}
assert allc == canon and len(allc) == 7, (allc, canon)
ald = lambda s: Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CX3H1](=O)'))
iod = lambda s: Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]'))
def alkane(s):   # C=O を CH2 に（還元→脱水→H2 付加で炭素骨格は保存）
    return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('(=O)', '').replace('=O', '')))
for k, s in carb.items():
    print(k, s, '非等価C', nC(s), 'アルデヒド', ald(s), 'ヨードホルム', iod(s), 'アルカン', alkane(s), nC(alkane(s)))
P = [k for k, s in carb.items() if nC(s) == 3 and ald(s)]; Q = [k for k, s in carb.items() if nC(s) == 3 and not ald(s)]
R = [k for k, s in carb.items() if nC(s) == 4 and iod(s)]
Tt = [k for k, s in carb.items() if nC(s) == 5 and not ald(s) and nC(alkane(s)) == 3]
print(P, Q, R, Tt)
assert P == ['2,2-ジメチルプロパナール'] and Q == ['3-ペンタノン'] and R == ['3-メチル-2-ブタノン'] and Tt == ['2-ペンタノン']
print('OK')
