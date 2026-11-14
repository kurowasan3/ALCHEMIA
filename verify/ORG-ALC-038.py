# ORG-ALC-038 イミン生成の平衡 検証
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import Counter
def c(s):
    return Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(s)).GetAtoms())
def M(s):
    A={'C':12,'H':1.0,'N':14,'O':16}; return sum(A[k]*v for k,v in c(s).items())
W=c('O')
# (1) アセトン + メチルアミン -> イミン + 水
assert c('CC(C)=O')+c('CN')==c('CC(C)=NC')+W; print(CalcMolFormula(Chem.MolFromSmiles('CC(C)=NC')))
# (4) 5-アミノペンタン-2-オン -> 2-メチル-1-ピロリン + 水
sub='CC(=O)CCCN'; prod='CC1=NCCC1'
assert c(sub)==c(prod)+W
print(CalcMolFormula(Chem.MolFromSmiles(sub)),M(sub),CalcMolFormula(Chem.MolFromSmiles(prod)),M(prod))
assert M(sub)==101 and M(prod)==83
mp=Chem.MolFromSmiles(prod); ri=mp.GetRingInfo(); assert ri.NumRings()==1 and ri.IsAtomInRingOfSize(1,5); print('ring size 5')
# 6-アミノヘキサン-2-オンなら6員環（比較用）
m6=Chem.MolFromSmiles('CC1=NCCCC1'); assert m6.GetRingInfo().IsAtomInRingOfSize(1,6)
# (5) 加水分解
assert c('C(=NCC)c1ccccc1')+W==c('O=Cc1ccccc1')+c('NCC')
# (2) イミンをつくれるアミン：N に H が2個
for name,s in [('メチルアミン','CN'),('ジメチルアミン','CNC'),('トリメチルアミン','CN(C)C'),('アニリン','Nc1ccccc1'),('アンモニア','N')]:
    m=Chem.AddHs(Chem.MolFromSmiles(s)); n=[a for a in m.GetAtoms() if a.GetSymbol()=='N'][0]
    print(name,'N-H 数',sum(1 for x in n.GetNeighbors() if x.GetSymbol()=='H'))
print('OK')
