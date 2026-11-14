# POL-PRT-071 トリプトファンからセロトニン・メラトニンへの変換 検証
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
from collections import Counter
A={'C':12,'H':1.0,'N':14,'O':16}
def c(s): return Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(s)).GetAtoms())
def M(s): return sum(A[k]*v for k,v in c(s).items())
S={'Trp':'NC(Cc1c[nH]c2ccccc12)C(=O)O','HTP':'NC(Cc1c[nH]c2ccc(O)cc12)C(=O)O',
   'Ser':'NCCc1c[nH]c2ccc(O)cc12','NAS':'CC(=O)NCCc1c[nH]c2ccc(O)cc12','Mel':'CC(=O)NCCc1c[nH]c2ccc(OC)cc12',
   'HIAL':'O=CCc1c[nH]c2ccc(O)cc12','HIAA':'OC(=O)Cc1c[nH]c2ccc(O)cc12'}
for k,s in S.items():
    m=Chem.MolFromSmiles(s); print(k,CalcMolFormula(m),M(s),FindMolChiralCenters(m,includeUnassigned=True))
mw={k:M(s) for k,s in S.items()}
assert [mw['HTP']-mw['Trp'],mw['Ser']-mw['HTP'],mw['NAS']-mw['Ser'],mw['Mel']-mw['NAS']]==[16,-44,42,14]
assert mw['Trp']==204 and mw['Mel']==232 and mw['HIAL']==175 and mw['HIAA']==191
# 元素分析 5-HIAA
x=c(S['HIAA']); m=mw['HIAA']
pc={k:round(A[k]*v/m*100,2) for k,v in x.items()}; print(pc)
assert abs(pc['C']-62.83)<0.01 and abs(pc['H']-4.71)<0.01 and abs(pc['N']-7.33)<0.01
# 組成式の決定
nC,nH,nN=62.83/12,4.71/1.0,7.33/14; nO=(100-62.83-4.71-7.33)/16
print([round(v/nN,2) for v in (nC,nH,nN,nO)])
# 収支：HTP -> Ser + CO2；Ser + CH3COOH -> NAS + H2O；Ser + O2 + H2O -> HIAL + NH3 + H2O2（参考）
assert c(S['HTP'])==c(S['Ser'])+c('O=C=O')
assert c(S['Ser'])+c('CC(=O)O')==c(S['NAS'])+c('O')
print('OK')
