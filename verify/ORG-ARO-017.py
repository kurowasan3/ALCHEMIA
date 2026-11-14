"""ORG-ARO-017 検証：C9H12O アルコールの性質"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, n_stereocenters, iodoform, has_cis_trans, can
c = {'ア': 'OCCCc1ccccc1', 'イ': 'OCC(C)c1ccccc1', 'ウ': 'CCC(O)c1ccccc1', 'エ': 'CC(O)Cc1ccccc1', 'オ': 'CC(C)(O)c1ccccc1'}
patt = Chem.MolFromSmarts('[CX4][OX2H1]')
cls = {}
for k, s in c.items():
    m = Chem.MolFromSmiles(s); assert formula(s) == 'C9H12O'
    ca = m.GetAtomWithIdx(m.GetSubstructMatch(patt)[0])
    cls[k] = ca.GetTotalNumHs()
print('carbinol H', cls)
assert [k for k in c if cls[k] == 2] == ['ア', 'イ'] and [k for k in c if cls[k] == 1] == ['ウ', 'エ'] and [k for k in c if cls[k] == 0] == ['オ']
st = [k for k in c if n_stereocenters(c[k])]; print('stereo', st); assert st == ['イ', 'ウ', 'エ']
io = [k for k in c if iodoform(c[k])]; print('iodoform', io); assert io == ['エ']
DEH = AllChem.ReactionFromSmarts('[C;!H0:1][C:2][OH:3]>>[C:1]=[C:2].[O:3]')
ct = []
for k, s in c.items():
    prods = set()
    for ps in DEH.RunReactants((Chem.MolFromSmiles(s),)):
        p = ps[0]
        try: Chem.SanitizeMol(p)
        except Exception: continue
        prods.add(Chem.MolToSmiles(p))
    print(k, prods, [has_cis_trans(p) for p in prods])
    if any(has_cis_trans(p) for p in prods): ct.append(k)
assert ct == ['ウ', 'エ']
P = 'C=C(C)c1ccccc1'
add = {'a': 'CC(C)c1ccccc1', 'b': 'BrCC(C)(Br)c1ccccc1', 'c': 'CC(C)(Br)c1ccccc1'}
r = {k: n_stereocenters(v) for k, v in add.items()}; print('adducts stereo', r)
assert r == {'a': 0, 'b': 1, 'c': 0}
report_equations('ORG-ARO-017')
