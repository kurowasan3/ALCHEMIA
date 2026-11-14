"""ORG-ARO-013 検証：一置換ベンゼンの C9H10O2 エステルの全列挙と A〜D の割り当て"""
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, can, has, iodoform, aldehyde, n_stereocenters
# 酸側は最後の原子が C(=O) に結合，アルコール側は最初の原子が O に結合
Racid = ['', 'C', 'CC', 'CC(C)', 'CCC', 'c1ccccc1', 'c1ccccc1C', 'c1ccccc1CC', 'c1ccccc1C(C)', 'CC(c1ccccc1)']
Ralc = ['C', 'CC', 'C(C)C', 'CCC', 'c1ccccc1', 'Cc1ccccc1', 'CCc1ccccc1', 'C(C)c1ccccc1']
esters = set()
for R, R2 in itertools.product(Racid, Ralc):
    smi = (f'{R}C(=O)O{R2}' if R else f'[H]C(=O)O{R2}')
    m = Chem.MolFromSmiles(smi)
    if m is None: continue
    s = Chem.MolToSmiles(m)
    if formula(s) == 'C9H10O2' and s.count('c') == 6 and len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts('c[!c;!#1]'))) == 1:
        esters.add(s)
print(len(esters), sorted(esters)); assert len(esters) == 6
HYD = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3][#6:4]>>[C:1](=[O:2])O.[O:3][#6:4]')
def hyd(s):
    ps = HYD.RunReactants((Chem.MolFromSmiles(s),))[0]
    out = []
    for p in ps: Chem.SanitizeMol(p); out.append(Chem.MolToSmiles(p))
    return out  # [acid, alcohol]
acidic = lambda s: has(s, 'C(=O)[OH]') or has(s, 'c[OH]')
cooh = lambda s: has(s, 'C(=O)[OH]')
phenol = lambda s: has(s, 'c[OH]')
benzoic_by_kmno4 = lambda s: has(s, 'c1ccccc1') and has(s, 'c[CX4;!H0]') and len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts('c[!c;!#1]'))) == 1
BENZOIC = can('OC(=O)c1ccccc1')
assign = {}
for e in esters:
    acid, alc = hyd(e)
    if cooh(acid) and aldehyde(acid) and not acidic(alc) and iodoform(alc): assign.setdefault('A', []).append(e)
    if phenol(alc) and cooh(acid): assign.setdefault('B', []).append(e)
    if cooh(acid) and not aldehyde(acid) and not has(acid, 'c') and not acidic(alc) and benzoic_by_kmno4(alc): assign.setdefault('C', []).append(e)
    if acid == BENZOIC and iodoform(alc): assign.setdefault('D', []).append(e)
print(assign)
exp = {'A': [can('O=COC(C)c1ccccc1')], 'B': [can('CCC(=O)Oc1ccccc1')], 'C': [can('CC(=O)OCc1ccccc1')], 'D': [can('CCOC(=O)c1ccccc1')]}
assert assign == exp
st = {k: n_stereocenters(v[0]) for k, v in assign.items()}; print('stereo', st)
assert st == {'A': 1, 'B': 0, 'C': 0, 'D': 0}
rest = sorted(esters - {v[0] for v in exp.values()}); print('others', rest, [hyd(r) for r in rest])
report_equations('ORG-ARO-013')
