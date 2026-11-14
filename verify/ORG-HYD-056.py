# ORG-HYD-056 検証：メチルヘプタンの不斉炭素，分子式，沸点データの内挿
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
mols={'octane':'CCCCCCCC','2-methylheptane':'CC(C)CCCCC','3-methylheptane':'CCC(C)CCCC','4-methylheptane':'CCCC(C)CCC','2,2,4-trimethylpentane':'CC(C)CC(C)(C)C'}
for k,s in mols.items():
    m=Chem.MolFromSmiles(s)
    cc=Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False)
    print(k,CalcMolFormula(m),cc)
    assert CalcMolFormula(m)=='C8H18'
assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles('CCC(C)CCCC'),includeUnassigned=True,useLegacyImplementation=False))==1
for s in ['CC(C)CCCCC','CCCC(C)CCC','CC(C)CC(C)(C)C','CCCCCCCC']:
    assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False))==0
# 直鎖アルカンの沸点（文献値，℃）：C6 68.7, C7 98.4, C8 125.7, C9 150.8, C10 174.1
bp={6:69,7:98,8:126,10:174}
est=(bp[8]+bp[10])/2
print('C9 推定（中点）',est,'実測 151')
assert abs(est-150.8)<2
# 増分が小さくなること
d=[bp[7]-bp[6],bp[8]-bp[7],(bp[10]-bp[8])/2]
print('増分',d); assert d[0]>d[1]>d[2]
print('OK')
