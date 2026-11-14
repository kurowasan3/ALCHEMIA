# POL-PRT-024 検証：SMILES から向かい合う 3 位置の D/A を機械的に判定し，塩基対の組合せと本数を全探索
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
B = {  # 名前: (SMILES（* = R）, 向かい合う位置の原子インデックス（ピリミジン 4,3,2 / プリン 6,1,2）, 種類)
 'C': ('Nc1ccn(*)c(=O)n1', (1, 8, 6), 'pyr'),
 'T': ('Cc1cn(*)c(=O)[nH]c1=O', (8, 7, 5), 'pyr'),
 'X': ('Nc1nc(=O)ccn1*', (3, 2, 1), 'pyr'),
 'G': ('Nc1nc2c(ncn2*)c(=O)[nH]1', (9, 11, 1), 'pur'),
 'A': ('Nc1ncnc2c1ncn2*', (1, 2, 3), 'pur'),
 'Z': ('Nc1nc(N)c2ncn(*)c2n1', (3, 2, 1), 'pur'),
 'Y': ('Nc1[nH]c(=O)nc2c1ncn2*', (1, 2, 3), 'pur'),
 "Y'": ('Nc1nc(O)nc2c1ncn2*', (1, 2, 3), 'pur'),
}
def role(m, idx):
    a = m.GetAtomWithIdx(idx)
    if a.GetSymbol() == 'N':          # 環の N
        return 'D' if a.GetTotalNumHs() > 0 else 'A'
    # 環の C：環外の置換基（*, 環原子以外）を見る
    subs = [n for n in a.GetNeighbors() if not n.IsInRing()]
    if not subs:
        return '-'                     # C-H
    s = subs[0]
    if s.GetSymbol() == 'N':
        return 'D'                     # -NH2
    if s.GetSymbol() == 'O':
        b = m.GetBondBetweenAtoms(a.GetIdx(), s.GetIdx())
        return 'A' if b.GetBondTypeAsDouble() == 2 else 'D'   # C=O は A，O-H は D
    return '?'
pat = {}
for k, (smi, idx, kind) in B.items():
    m = Chem.MolFromSmiles(smi)
    pat[k] = ''.join(role(m, i) for i in idx)
    print(k, pat[k], CalcMolFormula(Chem.MolFromSmiles(smi.replace('*', '[H]'))))
assert pat == {'C': 'DAA', 'T': 'ADA', 'X': 'AAD', 'G': 'ADD', 'A': 'DA-', 'Z': 'DAD', 'Y': 'DDA', "Y'": 'DAD'}
def pair(p, q):
    n = 0
    for a, b in zip(p, q):
        if a == '-' or b == '-':
            continue
        if a == b:
            return None
        n += 1
    return n
res = {(y, u): pair(pat[y], pat[u]) for y in 'CTX' for u in ['G', 'A', 'Z', 'Y']}
ok = {k: v for k, v in res.items() if v}
print(ok)
assert ok == {('C', 'G'): 3, ('T', 'A'): 2, ('T', 'Z'): 3, ('X', 'Y'): 3}
okY = {y: pair(pat[y], pat["Y'"]) for y in 'CTX'}
print("Y' :", okY); assert okY == {'C': None, 'T': 3, 'X': None}
# (2)
Zp = 30.0; G = (100 - 2*Zp)/2; assert G == 20.0
bp = 1.00e4; hz = bp*0.6*3 + bp*0.4*3; hn = bp*0.6*2 + bp*0.4*3
print(G, hz, hn); assert (hz, hn) == (3.0e4, 2.4e4)
print('OK')
