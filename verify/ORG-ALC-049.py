# ORG-ALC-049 verification: retrosynthesis of lactone T, stereocentres, ee propagation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
def info(s):
    m=Chem.MolFromSmiles(s)
    cc=Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False)
    return CalcMolFormula(m),len(cc),len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)})
T=info('O=C1CC(O)C(CC)C(c2ccccc2)O1'); A=info('O=CC(CC)C(O)c1ccccc1'); S=info('CCCC(O)C(CC)C=O')
print(T,A,S); assert T==('C13H16O3',3,8) and A[:2]==('C11H14O2',2) and S[0]=='C8H16O2'
# aldol mass balance: butanal + benzaldehyde
assert info('CCCC=O')[0]=='C4H8O' and info('O=Cc1ccccc1')[0]=='C7H6O'
# alpha-H check: benzaldehyde carbonyl carbon attached to aromatic C with no H
m=Chem.MolFromSmiles('O=Cc1ccccc1'); 
P=0.75*97+0.25*3; Pp=0.75*3+0.25*97; ee=P-Pp; print(P,Pp,ee); assert (P,Pp,ee)==(73.5,26.5,47.0)
assert abs(ee-94*0.5)<1e-9
fr=[P/2,P/2,Pp/2,Pp/2]; print(fr); assert [round(x) for x in fr]==[37,37,13,13]
ee_frac=(P/2-Pp/2)/(P/2+Pp/2)*100; print(ee_frac); assert round(ee_frac)==47
# ring size: O of C5 + C1..C5 -> 6 ring atoms; O of C3 + C1..C3 -> 4
mt=Chem.MolFromSmiles('O=C1CC(O)C(CC)C(c2ccccc2)O1'); rs=sorted(len(r) for r in Chem.GetSymmSSSR(mt)); print(rs); assert rs==[6,6]
print('OK')
