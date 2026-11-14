from collections import Counter
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def trees(n):
    if n==0: return [None]
    return [(a,b) for k in range(n) for a in trees(k) for b in trees(n-1-k)]   # (4-child, 6-child)
T=[len(trees(n)) for n in range(1,7)]; print(T); assert T==[1,2,5,14,42,132]
def chains(t):
    res=[]
    def walk(x):
        l=0
        while x is not None:
            l+=1
            if x[1] is not None: walk(x[1])
            x=x[0]
        res.append(l)
    walk(t); return tuple(sorted(res))
def rootchain(t):
    l=0
    while t is not None: l+=1; t=t[0]
    return l
def branched(t):
    return t is not None and ((t[0] is not None and t[1] is not None) or branched(t[0]) or branched(t[1]))
T5=trees(5)
assert sum(1 for t in T5 if not branched(t))==16
sel=[t for t in T5 if chains(t)==(1,2,2)]; assert len(sel)==10
rc=Counter(rootchain(t) for t in sel); print(rc); assert rc[1]==2 and rc[2]==8
# structures with root chain 1: root's 6-child is a maltose unit b-c; the other maltose d-e on b(6) or c(6)
r1=[t for t in sel if rootchain(t)==1]
exp={ (None, ((None,None),((None,None),None)) ), (None, ((None,((None,None),None)),None)) }
print(r1); assert set(r1)==exp
# build molecules with RDKit (alpha-D-Glc template) to make sure all 42 are distinct
def frag(c1,o4,o6): return f"{o6}C[C@H]1O[C@H]({c1})[C@H](O)[C@@H](O)[C@@H]1{o4}"
def build(t):
    frags=[]; ctr=[0]
    def rec(node, pm):
        ctr[0]+=1; me=ctr[0]
        c1='O' if pm is None else f'[*:{pm}]'; o4='O'; o6='O'
        if node[0] is not None: m=1000+me; o4=f'O[*:{m}]'; rec(node[0],m)
        if node[1] is not None: m=2000+me; o6=f'[*:{m}]O'; rec(node[1],m)
        frags.append(frag(c1,o4,o6))
    rec(t,None)
    return Chem.MolToSmiles(Chem.molzip(Chem.MolFromSmiles('.'.join(frags))))
assert build((None,None))==Chem.CanonSmiles('C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O')
assert build(((None,None),None))==Chem.CanonSmiles('C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O')
S=[build(t) for t in T5]; assert len(set(S))==42
assert {CalcMolFormula(Chem.MolFromSmiles(s)) for s in S}=={'C30H52O26'}
print('POL-SUG-028 OK')
