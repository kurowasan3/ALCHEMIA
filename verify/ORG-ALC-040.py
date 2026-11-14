# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/out/items/ORG-ALC-040/item.json'
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
# 簡単のため手で列挙（ブテン骨格の各 H を Ph に置換した構造，重複なし）
cands=['C=CCCc1ccccc1','C=CC(C)c1ccccc1','C=C(CC)c1ccccc1','C(=CCC)c1ccccc1',
       'CC=CCc1ccccc1','CC=C(C)c1ccccc1','C=C(C)Cc1ccccc1','CC(C)=Cc1ccccc1']
rxn=AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.[C:2]=O')
hits=[]
for s in cands:
    m=Chem.MolFromSmiles(s); assert CalcMolFormula(m)=='C10H12', s
    ps=rxn.RunReactants((m,))
    for p in ps[0]: Chem.SanitizeMol(p)
    prods=sorted(Chem.MolToSmiles(p) for p in ps[0])
    if 'C=O' in prods:
        other=[p for p in prods if p!='C=O'][0]
        km=Chem.MolFromSmiles(other)
        ketone=km.HasSubstructMatch(Chem.MolFromSmarts('[#6][CX3](=O)[#6]'))
        iodo=km.HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX3](=O)'))
        print(s, prods, 'ketone',ketone,'iodoform',iodo)
        if ketone: hits.append((s,other,iodo))
assert len(hits)==2
assert sum(h[2] for h in hits)==1
A=[h for h in hits if h[2]][0]; Ap=[h for h in hits if not h[2]][0]
assert Chem.CanonSmiles(A[0])==Chem.CanonSmiles('C=C(C)Cc1ccccc1')
assert Chem.CanonSmiles(Ap[0])==Chem.CanonSmiles('C=C(CC)c1ccccc1')
def nchir(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False))
assert nchir('CC(O)Cc1ccccc1')==1 and nchir('CC(C)Cc1ccccc1')==0 and nchir('CCC(C)c1ccccc1')==1
print('不斉炭素 OK')

print('ALL OK')
