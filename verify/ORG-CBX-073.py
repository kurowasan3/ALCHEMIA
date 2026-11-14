# ORG-CBX-073 ベンジルエステル・tert-ブチルエステルの選択的切断 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def mw(smi):
    m=Chem.MolFromSmiles(smi); f=CalcMolFormula(m)
    A={'C':12,'H':1.0,'O':16}
    mass=sum(A[a.GetSymbol()] for a in Chem.AddHs(m).GetAtoms())
    return f,mass
S={'P':'O=C(OCc1ccccc1)/C=C/C(=O)OC(C)(C)C','monoBn_fum':'O=C(O)/C=C/C(=O)OCc1ccccc1',
   'monotBu_succ':'OC(=O)CCC(=O)OC(C)(C)C','fumaric':'OC(=O)/C=C/C(=O)O','toluene':'Cc1ccccc1',
   'isobutene':'C=C(C)C','BnOH':'OCc1ccccc1','tBuOH':'CC(C)(C)O','monotBu_fum':'OC(=O)/C=C/C(=O)OC(C)(C)C'}
R={k:mw(v) for k,v in S.items()}
for k,v in R.items(): print(k,v)
assert R['P']==('C15H18O4',262.0) and R['monoBn_fum']==('C11H10O4',206.0)
# 原子の収支（Y: P -> monoBn + isobutene ; X/H2: P + 2H2 -> monotBu_succ + toluene）
from collections import Counter
def cnt(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); return Counter(a.GetSymbol() for a in m.GetAtoms())
assert cnt(S['P'])==cnt(S['monoBn_fum'])+cnt(S['isobutene'])
assert cnt(S['P'])+Counter({'H':4})==cnt(S['monotBu_succ'])+cnt(S['toluene'])
assert cnt(S['P'])+Counter({'H':4,'O':2})==cnt(S['fumaric'])+cnt(S['BnOH'])+cnt(S['tBuOH'])
# 2.62 g の P から
n=2.62/262; print(n*206); assert round(n*206,2)==2.06
print('OK')
