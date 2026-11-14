# POL-SUG-046 検証：単糖の鎖状構造と環状構造の立体異性体の数
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def n(smi):
    m=Chem.MolFromSmiles(smi)
    opts=StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    iso={Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m,options=opts)}
    return len(iso), len(Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False))
S={'aldohexose chain':'OCC(O)C(O)C(O)C(O)C=O','aldohexose pyranose':'OCC1OC(O)C(O)C(O)C1O','aldohexose furanose':'OCC(O)C1OC(O)C(O)C1O',
   'aldopentose chain':'OCC(O)C(O)C(O)C=O','aldopentose pyranose':'OC1COC(O)C(O)C1O','aldopentose furanose':'OCC1OC(O)C(O)C1O',
   'ketohexose chain':'OCC(=O)C(O)C(O)C(O)CO','ketohexose furanose':'OCC1(O)OC(CO)C(O)C1O','ketohexose pyranose':'OCC1(O)OCC(O)C(O)C1O'}
R={}
for k,v in S.items():
    R[k]=n(v); print(f'{k:22s}', CalcMolFormula(Chem.MolFromSmiles(v)), 'stereoisomers,centers', R[k])
assert R['aldohexose chain']==(16,4) and R['aldohexose pyranose']==(32,5) and R['aldohexose furanose']==(32,5)
assert R['aldopentose chain']==(8,3) and R['aldopentose pyranose']==(16,4) and R['aldopentose furanose']==(16,4)
assert R['ketohexose chain']==(8,3) and R['ketohexose furanose']==(16,4)
# α形とβ形は5つの不斉中心のうち1位だけが逆→鏡像異性体ではなくジアステレオマー
print('alpha/beta: 1 of 5 centers inverted -> diastereomers')
