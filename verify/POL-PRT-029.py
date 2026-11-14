from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
from collections import Counter
import re
def cnt(f): 
    c=Counter()
    for el,n in re.findall(r'([A-Z][a-z]?)(\d*)',f): c[el]+=int(n or 1)
    return c
def fmt(c): return ''.join(f'{e}{c[e] if c[e]>1 else ""}' for e in ['C','H','N','O','P'] if c[e])
smi={'adenine':'Nc1ncnc2[nH]cnc12','thymine':'Cc1c[nH]c(=O)[nH]c1=O','dRib':'OCC1OC(O)CC1O','Rib':'OCC1OC(O)C(O)C1O',
     'dAMP':'Nc1ncnc2n(cnc12)C1CC(O)C(COP(=O)(O)O)O1','dTMP':'Cc1cn(C2CC(O)C(COP(=O)(O)O)O2)c(=O)[nH]c1=O',
     'pdApdT':'Nc1ncnc2n(cnc12)C1CC(OP(=O)(O)OCC3OC(n4cc(C)c(=O)[nH]c4=O)CC3O)C(COP(=O)(O)O)O1'}
f={k:F(Chem.MolFromSmiles(v)) for k,v in smi.items()}; print(f)
assert f['Rib']=='C5H10O5' and f['dRib']=='C5H10O4' and f['adenine']=='C5H5N5' and f['thymine']=='C5H6N2O2'
W=cnt('H2O'); Pa=cnt('H3PO4')
dAMP=cnt(f['adenine'])+cnt(f['dRib'])+Pa; dAMP.subtract(W); dAMP.subtract(W)
assert fmt(dAMP)==f['dAMP']=='C10H14N5O6P'
dTMP=cnt(f['thymine'])+cnt(f['dRib'])+Pa; dTMP.subtract(W); dTMP.subtract(W)
assert fmt(dTMP)==f['dTMP']=='C10H15N2O8P'
di=dAMP+dTMP; di.subtract(W)
assert fmt(di)==f['pdApdT']=='C20H27N7O13P2'
n=10; assert 2*n+(n-1)==29==3*n-1
print('POL-PRT-029 OK')
