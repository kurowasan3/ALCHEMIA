"""ORG-HYD-016 検証：C5H8 異性体の判別"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, n_stereocenters, can, has, iodoform
c = {'a': 'C#CCCC', 'b': 'CC#CCC', 'c': 'C#CC(C)C', 'd': 'C=CC=CC', 'e': 'C=C(C)C=C', 'f': 'C=CCC=C'}
br2 = {'a': 'BrC(Br)C(Br)(Br)CCC', 'b': 'CC(Br)(Br)C(Br)(Br)CC', 'c': 'BrC(Br)C(Br)(Br)C(C)C',
       'd': 'BrCC(Br)C(Br)C(Br)C', 'e': 'BrCC(C)(Br)C(Br)CBr', 'f': 'BrCC(Br)CC(Br)CBr'}
def hydro(s): return can(s.replace('#', '').replace('=', ''))
term = {k: has(s, '[CH]#C') for k, s in c.items()}
for k, s in c.items():
    assert formula(s) == 'C5H8'
    assert formula(br2[k]) == 'C5H8Br4' and hydro(br2[k].replace('(Br)', '').replace('Br', '')) == hydro(s)
    print(k, s, 'terminal', term[k], 'H2->', hydro(s), 'Br4 stereo', n_stereocenters(br2[k]))
assert [k for k in c if term[k]] == ['a', 'c']
ib = can('CC(C)CC')
P = [k for k in c if term[k] and hydro(c[k]) == ib]; Q = [k for k in c if not term[k] and hydro(c[k]) == ib]
R = [k for k in c if not term[k] and n_stereocenters(br2[k]) == 0]
print('P', P, 'Q', Q, 'R', R); assert (P, Q, R) == (['c'], ['e'], ['b'])
ket = 'CC(=O)C(C)C'; assert formula(ket) == 'C5H10O' and iodoform(ket)
r12 = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](Br)[C:2]Br')
r14 = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[C:3]=[C:4]>>[C:1](Br)-[C:2]=[C:3]-[C:4]Br')
prods = set()
m = Chem.MolFromSmiles(c['e'])
for rx in (r12, r14):
    for ps in rx.RunReactants((m,)):
        p = ps[0]; Chem.SanitizeMol(p); prods.add(Chem.MolToSmiles(p))
print('C5H8Br2 from isoprene', prods); assert len(prods) == 3
bd = set()
for rx in (r12, r14):
    for ps in rx.RunReactants((Chem.MolFromSmiles('C=CC=C'),)):
        p = ps[0]; Chem.SanitizeMol(p); bd.add(Chem.MolToSmiles(p))
print('butadiene comparison', bd); assert len(bd) == 2
assert n_stereocenters(br2['d']) == 3
report_equations('ORG-HYD-016')
