# ORG-STR-071 検証：分子量 44（H=1,C=12,N=14,O=16）の分子式の列挙，価標の偶奇，C2H4O の構造異性体
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
cands=[]
for a in range(0,5):
    for c in range(0,5):
        for d in range(0,4):
            b=44-12*a-14*c-16*d
            if b<0: continue
            if a+b+c+d>12 or a+b+c+d<2: continue
            cands.append((a,b,c,d))
print('候補', cands)
assert len(cands)==8
odd=[x for x in cands if (4*x[0]+x[1]+3*x[2]+2*x[3])%2==1]
print('価標の総数が奇数', odd)
assert sorted(odd)==sorted([(0,2,3,0),(1,2,1,1),(2,6,1,0)])
# 偶奇は b+c の偶奇と一致
assert all(((4*a+b+3*c+2*d)%2)==((b+c)%2) for a,b,c,d in cands)
# 不飽和度 >= 0
for a,b,c,d in cands:
    U=(2*a-b+c)/2+1
    assert U>=0, (a,b,c,d)
# C2H4O の構造（環・鎖すべて）を RDKit で確認
smis=['CC=O','C=CO','C1CO1','C1OC1']
can=set(Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smis)
print('C2H4O の構造', can)
assert len(can)==3
for s in can:
    assert CalcMolFormula(Chem.MolFromSmiles(s))=='C2H4O'
# 2個の重原子 C,C と O の結合のしかたを網羅：O の位置（末端 C-C-O, C-O-C, 環）と二重結合の位置で 3 通り＋カルベン等を除く
print('OK')
