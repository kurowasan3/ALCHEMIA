# ORG-STR-056 検証：ピリジン環の置換体の数（ベンゼンとの比較）と隣接ジカルボン酸
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import permutations, combinations
def build(ring, subs):
    # ring: 'n' or 'c' at position 0, subs: dict pos->smiles
    atoms=[]
    for i in range(6):
        a=ring if i==0 else 'c'
        if i==0: a+='1'
        if i in subs: a+='(%s)'%subs[i]
        atoms.append(a)
    return Chem.CanonSmiles(''.join(atoms)+'1')
def count(ring, groups):
    pos=range(1,6) if ring=='n' else range(6)
    out=set()
    for ps in permutations(pos,len(groups)):
        out.add(build(ring, dict(zip(ps,groups))))
    return out
res={}
res['py ClMe']=count('n',['Cl','C']); res['py Cl2']=count('n',['Cl','Cl']); res['py Me']=count('n',['C'])
res['bz ClMe']=count('c',['Cl','C']); res['bz Cl2']=count('c',['Cl','Cl'])
res['py COOH2']=count('n',['C(=O)O','C(=O)O'])
for k,v in res.items(): print(k, len(v), {CalcMolFormula(Chem.MolFromSmiles(s)) for s in v})
assert [len(res[k]) for k in ['py ClMe','py Cl2','py Me','bz ClMe','bz Cl2','py COOH2']]==[10,6,3,3,3,6]
# 隣り合う COOH（環の隣接炭素）をもつもの
patt=Chem.MolFromSmarts('c(C(=O)O):c(C(=O)O)')
adj=[s for s in res['py COOH2'] if Chem.MolFromSmiles(s).HasSubstructMatch(patt)]
print('adjacent', adj); assert len(adj)==2
