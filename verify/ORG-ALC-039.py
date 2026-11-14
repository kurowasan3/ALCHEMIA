# ORG-ALC-039 環状ヘミアセタールと 1,2-プロパンジオールの混合アセタール
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
S={'5-ヒドロキシペンタナール':'OCCCCC=O','環状ヘミアセタール':'OC1CCCCO1','1,2-プロパンジオール':'CC(O)CO',
 'X':'CC(O)COC1CCCCO1','Y':'OCC(C)OC1CCCCO1','Z':'OCCCCC1OCC(C)O1'}
iodo=Chem.MolFromSmarts('[CH3][CX4;H1]([OX2H1])[#6,#1]')
acetal=Chem.MolFromSmarts('[CX4;H1](O[#6])(O[#6])')
for k,s in S.items():
    m=Chem.MolFromSmiles(s); ri=m.GetRingInfo()
    print(k,CalcMolFormula(m),'rings',[len(r) for r in ri.AtomRings()],'iodoform',m.HasSubstructMatch(iodo),'acetal C',m.HasSubstructMatch(acetal),
          'free OH',len(m.GetSubstructMatches(Chem.MolFromSmarts('[OX2H1]'))))
# 加水分解生成物が同じ：アセタール C を切ると 5-ヒドロキシペンタナール + プロパンジオール
print('C5H10O2 + C3H8O2 - H2O =', 'C8H16O3')
