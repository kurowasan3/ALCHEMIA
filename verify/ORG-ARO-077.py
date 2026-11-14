# ORG-ARO-077 検証：OH・OCH3・CH3 の三置換ベンゼンの異性体と，メチル化による収束
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import permutations
def ring(subs):  # subs: list of 6 substituent strings or None (positions 1..6)
    s=''
    atoms=[]
    for i,x in enumerate(subs):
        a='c' if i else 'c1'
        if x: a+='(%s)'%x
        atoms.append(a)
    return Chem.CanonSmiles(''.join(atoms)+'1') if False else Chem.CanonSmiles(atoms[0]+''.join(atoms[1:])+'1')
phen=set()
for pos in permutations(range(6),3):
    subs=[None]*6
    subs[pos[0]]='O'; subs[pos[1]]='OC'; subs[pos[2]]='C'
    phen.add(ring(subs))
print('phenols', len(phen), {CalcMolFormula(Chem.MolFromSmiles(s)) for s in phen})
assert len(phen)==10
# メチル化：OH -> OCH3
rxn_map={}
for s in phen:
    m=Chem.MolFromSmiles(s)
    p=Chem.ReplaceSubstructs(m, Chem.MolFromSmarts('[OX2H]'), Chem.MolFromSmiles('OC'), replaceAll=True)[0]
    ps=Chem.CanonSmiles(Chem.MolToSmiles(p))
    rxn_map.setdefault(ps,[]).append(s)
print('methylated products', len(rxn_map), {CalcMolFormula(Chem.MolFromSmiles(s)) for s in rxn_map})
assert len(rxn_map)==6
for k,v in sorted(rxn_map.items(), key=lambda kv: len(kv[1])):
    print(len(v), k, v)
# 酸化（CH3->COOH）→脱炭酸（COOH->H）で 1,3-ジメトキシベンゼンになるもの
target=Chem.CanonSmiles('COc1cccc(OC)c1')
hits=[]
for k,v in rxn_map.items():
    m=Chem.MolFromSmiles(k)
    p=Chem.ReplaceSubstructs(m, Chem.MolFromSmarts('[CH3;!$(CO)]'), Chem.MolFromSmiles('[H]'), replaceAll=True)[0]
    ps=Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(p)))
    if ps==target: hits.append((k,len(v)))
print('to m-dimethoxybenzene:', hits)
two=[h for h in hits if h[1]==2]
assert len(two)==1
R=two[0][0]; print('R =',R, CalcMolFormula(Chem.MolFromSmiles(R)))
print('P,Q =', rxn_map[R])
# S: R の CH3 -> COOH
S=Chem.ReplaceSubstructs(Chem.MolFromSmiles(R), Chem.MolFromSmarts('[CH3;!$(CO)]'), Chem.MolFromSmiles('C(=O)O'))[0]
print('S', CalcMolFormula(Chem.MolFromSmiles(Chem.MolToSmiles(S))))
assert Chem.CanonSmiles(R)==Chem.CanonSmiles('Cc1ccc(OC)cc1OC')
assert set(rxn_map[R])=={Chem.CanonSmiles('Cc1ccc(OC)cc1O'),Chem.CanonSmiles('Cc1ccc(O)cc1OC')}
