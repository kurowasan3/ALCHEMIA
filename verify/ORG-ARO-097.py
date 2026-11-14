# ORG-ARO-097 verification: model of dehydroaromatization of alkanes (6 consecutive carbons close a benzene ring)
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx

def products(smi):
    m=Chem.MolFromSmiles(smi)
    G=nx.Graph(); G.add_edges_from((b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds())
    out=set()
    for s in G.nodes:
        for t in G.nodes:
            if s<t:
                for p in nx.all_simple_paths(G,s,t,cutoff=5):
                    if len(p)!=6: continue
                    # ring closure s-t; each ring carbon must end up with <=3 carbon neighbours
                    if any(G.degree(a)+(1 if a in (s,t) else 0)>3 for a in p): continue
                    rw=Chem.RWMol(m)
                    rw.AddBond(s,t,Chem.BondType.SINGLE)
                    ring=p
                    for i in range(6):
                        a,b=ring[i],ring[(i+1)%6]
                        rw.GetBondBetweenAtoms(a,b).SetBondType(Chem.BondType.DOUBLE if i%2==0 else Chem.BondType.SINGLE)
                    mm=rw.GetMol(); Chem.SanitizeMol(mm)
                    out.add(Chem.MolToSmiles(mm))
    return out
for name,smi in [('octane','CCCCCCCC'),('decane','CCCCCCCCCC'),('undecane','CCCCCCCCCCC'),('tridecane','CCCCCCCCCCCCC'),
                 ('2-methylheptane','CC(C)CCCCC'),('2,2-dimethylbutane','CC(C)(C)CC'),('2,2-dimethylhexane','CC(C)(C)CCCC'),('3-methylheptane','CCC(C)CCCC')]:
    ps=products(smi)
    print(name,len(ps),[ (p,CalcMolFormula(Chem.MolFromSmiles(p))) for p in sorted(ps)])
P=products('CCCCCCCCCCC'); assert len(P)==3
assert Chem.CanonSmiles('CCCCCc1ccccc1') in P and Chem.CanonSmiles('CCCCc1ccccc1C') in P and Chem.CanonSmiles('CCCc1ccccc1CC') in P
P13=products('CCCCCCCCCCCCC'); assert len(P13)==4 and Chem.CanonSmiles('CCCCc1ccccc1CCC') in P13
assert products('CC(C)CCCCC')=={Chem.CanonSmiles('Cc1cccc(C)c1')}
assert products('CC(C)(C)CC')==set() and products('CC(C)(C)CCCC')==set()
assert check_equation('C11H24 + 4C6H12 -> C11H16 + 4C6H14') is None
assert check_equation('C11H24 -> C11H16 + 4H2') is None
# mass of acceptor per 1.00 g undecane: 4*84/156
print('acceptor g per g undecane',4*84.0/156.0)
assert round(4*84.0/156.0,2)==2.15
print('OK')
