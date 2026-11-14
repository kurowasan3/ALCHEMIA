# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/outh/items/ORG-HYD-055/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
X='C1CCC2=C(C1)CCCC2'; Xp='C1CCC2=CCCCC2C1'
for s,f in [(X,'C10H16'),(Xp,'C10H16')]: assert CalcMolFormula(Chem.MolFromSmiles(s))==f
rxn=AllChem.ReactionFromSmarts('[C:1]=[C:2]>>([C:1]=O.[C:2]=O)')
def ozo(s):
    ps=rxn.RunReactants((Chem.MolFromSmiles(s),))
    out=set()
    for p in ps:
        for q in p: Chem.SanitizeMol(q); out.add(Chem.MolToSmiles(q))
    return out
Y=ozo(X); Z=ozo(Xp); print(Y,Z)
assert Y=={Chem.CanonSmiles('O=C1CCCCC(=O)CCCC1')}
assert Z=={Chem.CanonSmiles('O=CCCCC1CCCCC1=O')}
for s in list(Y)+list(Z): assert CalcMolFormula(Chem.MolFromSmiles(s))=='C10H16O2'
def nCl(k):
    ring=['C']*10; ring[0]='C(=O)'; ring[k-1]='C(=O)'
    smi='C1'+''.join(ring[1:])+'1'; smi='O=C1'+''.join('C(=O)' if i==k-1 else 'C' for i in range(1,10))+'1'
    m=Chem.MolFromSmiles(smi); assert CalcMolFormula(m)=='C10H16O2'
    ranks=list(Chem.CanonicalRankAtoms(m,breakTies=False))
    cls={ranks[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol()=='C' and a.GetTotalNumHs()>0}
    return len(cls)
res={k:nCl(k) for k in (3,4,5,6)}; print('一塩素置換体', res)
assert res=={3:5,4:4,5:5,6:2}
ch=lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False))
assert ch(list(Z)[0])==1 and ch(list(Y)[0])==0 and ch(X)==0
svg=open(ITEM.replace('item.json','solution.svg')).read(); assert 'viewBox' in svg and '<script' not in svg

print('ALL OK')
