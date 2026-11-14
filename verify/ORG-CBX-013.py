"""ORG-CBX-013 検証：ヒドロキシ酸の不斉炭素・還元体・ラクトンの環の員数"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from g09lib import report_equations, formula, n_stereocenters
acids = {'ア': 'OCCCC(=O)O', 'イ': 'CC(O)CC(=O)O', 'ウ': 'CCC(O)C(=O)O', 'エ': 'OCC(C)C(=O)O', 'オ': 'CC(C)(O)C(=O)O'}
red = {k: v.replace('C(=O)O', 'CO') for k, v in acids.items()}
patt = Chem.MolFromSmarts('[CX4][OX2H1]')
cpat = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
res = {}
for k, s in acids.items():
    m = Chem.MolFromSmiles(s)
    assert formula(s) == 'C4H8O3'
    o = m.GetSubstructMatch(patt)[1]
    c = m.GetSubstructMatch(cpat)[0]
    ring = len(Chem.GetShortestPath(m, o, c))
    res[k] = (n_stereocenters(s), n_stereocenters(red[k]), ring)
    print(k, s, 'stereo', res[k][0], 'reduced', red[k], res[k][1], 'ring', ring)
assert [k for k in res if res[k][0]] == ['イ', 'ウ', 'エ']
assert [k for k in res if res[k][0] and not res[k][1]] == ['エ']
assert [k for k in res if res[k][2] == 5] == ['ア']
gbl = 'O=C1CCCO1'; lactide = 'CC1OC(=O)C(C)OC1=O'
print('GBL', formula(gbl), 'lactide', formula(lactide))
assert formula(gbl) == 'C4H6O2' and formula(lactide) == 'C6H8O4'
mlac = Chem.MolFromSmiles(lactide); ri = mlac.GetRingInfo()
print('lactide ring sizes', [len(r) for r in ri.AtomRings()]); assert [len(r) for r in ri.AtomRings()] == [6]
report_equations('ORG-CBX-013')
