# ORG-ARO-078 検証：アニスアルデヒドの2つの合成経路
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import AllChem
M={'H':1.0,'C':12.0,'O':16.0}
def mw(s):
    m=Chem.AddHs(Chem.MolFromSmiles(s)); return sum(M[a.GetSymbol()] for a in m.GetAtoms())
S={'phenol':'Oc1ccccc1','A':'OCc1ccc(O)cc1','B':'O=Cc1ccc(O)cc1','anisaldehyde':'COc1ccc(C=O)cc1',
   'trans-anethole':'COc1ccc(/C=C/C)cc1','cis-anethole':'COc1ccc(/C=C\\C)cc1','estragole':'COc1ccc(CC=C)cc1',
   'acetaldehyde':'CC=O','homoanisaldehyde':'COc1ccc(CC=O)cc1','HCHO':'C=O'}
for k,v in S.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(v)), mw(v))
# オゾン分解（C=C を切って両側を C=O に）
rxn=AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.[C:2]=O')
def ozon(s):
    m=Chem.MolFromSmiles(s)
    ps=rxn.RunReactants((m,))
    out=set()
    for p in ps:
        frags=[]
        for x in p:
            Chem.SanitizeMol(x); frags.append(Chem.CanonSmiles(Chem.MolToSmiles(x,isomericSmiles=False)))
        out.add(tuple(sorted(frags)))
    return out
# ベンゼン環の芳香族 C=C は SMARTS の '=' (非芳香族二重結合) に一致しない
for k in ['trans-anethole','cis-anethole','estragole']:
    print(k,'->',ozon(S[k]))
assert ozon(S['trans-anethole'])==ozon(S['cis-anethole'])=={tuple(sorted([Chem.CanonSmiles(S['anisaldehyde']),Chem.CanonSmiles('CC=O')]))}
assert ozon(S['estragole'])=={tuple(sorted([Chem.CanonSmiles(S['homoanisaldehyde']),Chem.CanonSmiles('C=O')]))}
# エストラゴールにシス-トランス異性体がない
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
print('estragole stereo', len(list(EnumerateStereoisomers(Chem.MolFromSmiles('COc1ccc(CC=C)cc1')))), 'anethole stereo', len(list(EnumerateStereoisomers(Chem.MolFromSmiles('COc1ccc(C=CC)cc1')))))
m=1.0*mw(S['anisaldehyde'])/mw(S['trans-anethole']); print('mass kg', m); assert round(m,2)==0.92
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
print(check_equation('C6H5OH + HCHO -> C6H4(OH)CH2OH'))
