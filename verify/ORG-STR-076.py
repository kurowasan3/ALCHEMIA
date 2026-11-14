# ORG-STR-076 verification: CIP labels change on hydrogenation without touching the stereocentre; SN2 inversion
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
def lab(smi):
    m=Chem.MolFromSmiles(smi); rdCIPLabeler.AssignCIPLabels(m)
    return [a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode')]
R_ene='C=C[C@H](C)CCC'; assert lab(R_ene)==['R']           # (R)-3-methyl-1-hexene
hyd=R_ene.replace('C=C','CC',1); print('hydrogenated',hyd,lab(hyd)); assert lab(hyd)==['S']
# product is chiral (ethyl vs propyl): enantiomers differ in canonical isomeric smiles
m1=Chem.MolToSmiles(Chem.MolFromSmiles(hyd)); m2=Chem.MolToSmiles(Chem.MolFromSmiles('CC[C@@H](C)CCC')); assert m1!=m2
R_br='CC[C@H](Br)C'; assert lab(R_br)==['R']
ret=R_br.replace('Br','O'); assert lab(ret)==['R']          # retention would keep R
inv='CC[C@@H](O)C'; print('inversion',lab(inv)); assert lab(inv)==['S']
for e in ['CH2=CHCH(CH3)CH2CH2CH3 + H2 -> CH3CH2CH(CH3)CH2CH2CH3','CH3CH2CHBrCH3 + OH- -> CH3CH2CH(OH)CH3 + Br-']:
    r=check_equation(e); print(e,'->',r or 'OK'); assert r is None
print('OK')
