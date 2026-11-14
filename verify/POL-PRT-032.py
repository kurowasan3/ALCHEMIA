import itertools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
smi={'Gly':'NCC(=O)O','Cys':'SCC(N)C(=O)O','Phe':'c1ccccc1CC(N)C(=O)O','Leu':'CC(C)CC(N)C(=O)O','Ile':'CCC(C)C(N)C(=O)O'}
CH={k:len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(v),includeUnassigned=True,useLegacyImplementation=False)) for k,v in smi.items()}
print(CH); assert CH=={'Gly':0,'Cys':1,'Phe':1,'Leu':1,'Ile':2}
assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi['Leu']))=='C6H13NO2'==rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi['Ile']))
def ch(f): return sum(CH[a] for a in f)
def sul(f): return 'Cys' in f
def frags(s): return {L:[s[:L],s[5-L:]] for L in range(1,5)}
def obs(s, use_s):
    o=[]
    for L,fs in frags(s).items():
        o.append((L,tuple(sorted((ch(f),sul(f)) if use_s else (ch(f),) for f in fs))))
    return tuple(o)
target=('Cys','Gly','Ile','Phe','Gly')
assert ch(target)==4
# table data from the problem
table={'A':(4,True),'B':(3,False),'C':(3,True),'D':(3,False),'E':(1,True),'F':(1,False),'G':(0,False),'H':(1,True)}
fr=frags(target)
assert sorted((ch(f),sul(f)) for f in fr[4])==sorted([table['A'],table['B']])
assert sorted((ch(f),sul(f)) for f in fr[3])==sorted([table['C'],table['D']])
assert sorted((ch(f),sul(f)) for f in fr[2])==sorted([table['E'],table['F']])
assert sorted((ch(f),sul(f)) for f in fr[1])==sorted([table['G'],table['H']])
allf=[f for L in fr for f in fr[L]]; assert len(set(allf))==8
def cands(use_s, nterm=None):
    c=set()
    for x in ('Leu','Ile'):
        if 2+CH[x]!=4: continue
        for p in set(itertools.permutations(['Gly','Gly','Cys','Phe',x])):
            if obs(p,use_s)==obs(target,use_s) and (nterm is None or p[0]==nterm): c.add(p)
    return c
# Leu excluded by total chirality
assert 2+CH['Leu']!=4
c1=cands(False); print(len(c1),c1); assert len(c1)==4
c2=cands(True); print(len(c2),c2); assert len(c2)==2
c3=cands(True,'Cys'); print(c3); assert c3=={target}
# pairing: A(4) with G(0), B(3) with H(1)
assert table['A'][0]+table['G'][0]==4 and table['B'][0]+table['H'][0]==4
print('POL-PRT-032 OK')
