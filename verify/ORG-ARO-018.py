"""ORG-ARO-018 検証：C8H8O 異性体と側鎖酸化"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from g09lib import report_equations, formula, iodoform, aldehyde, has, can
S = dict(P='Cc1ccc(C=O)cc1', Q='O=CCc1ccccc1', R='CC(=O)c1ccccc1', S='C=Cc1ccccc1O', T='CCc1ccccc1O')
for k, s in S.items(): print(k, formula(s))
assert all(formula(S[k]) == 'C8H8O' for k in 'PQRS') and formula(S['T']) == 'C8H10O'
def side_oxidize(smi):
    """ベンゼン環に直接結合した炭素側鎖（H をもつ炭素，またはアルデヒド）をすべて COOH に"""
    m = Chem.MolFromSmiles(smi)
    ring = [a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()]
    subs = []
    for i in ring:
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            if not nb.GetIsAromatic() and nb.GetSymbol() == 'C':
                subs.append(i)
    rw = Chem.RWMol(Chem.MolFromSmiles('c1ccccc1'))
    # 置換位置だけ取り出して COOH を付けた分子を作る
    pos = sorted(ring.index(i) for i in subs)
    atoms = ['c'] * 6
    for p in pos: atoms[p] = 'c(C(=O)O)'
    return can(atoms[0] + '1' + ''.join(atoms[1:]) + '1')
print('P ox', side_oxidize(S['P']), 'Q ox', side_oxidize(S['Q']))
assert side_oxidize(S['P']) == can('OC(=O)c1ccc(C(=O)O)cc1') and side_oxidize(S['Q']) == can('OC(=O)c1ccccc1')
assert aldehyde(S['P']) and aldehyde(S['Q']) and not aldehyde(S['R']) and iodoform(S['R'])
assert has(S['S'], 'c[OH]') and can(S['S'].replace('C=C', 'CC')) == can(S['T'])  # S + H2 -> T（o-エチルフェノール）
xyl = {'o': (-25.2, 144.4), 'm': (-47.9, 139.1), 'p': (13.3, 138.4)}
assert max(xyl, key=lambda k: xyl[k][0]) == 'p' and max(xyl, key=lambda k: xyl[k][1]) == 'o'
report_equations('ORG-ARO-018')
