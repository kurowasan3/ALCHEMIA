"""ORG-STR-054 検証：シス/トランスの判定・キラリティー・モデルのエネルギー・存在比"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def cis_trans(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 7; ps.enforceChirality = True
    assert AllChem.EmbedMolecule(m, ps) == 0
    Chem.AssignStereochemistryFrom3D(m)
    assert Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(m))) == Chem.CanonSmiles(smi), 'stereo lost'
    conf = m.GetConformer(); ri = m.GetRingInfo().AtomRings()[0]
    P = np.array([list(conf.GetAtomPosition(i)) for i in ri]); c = P.mean(0)
    n = np.linalg.svd(P - c)[2][2]
    pos = lambda i: np.array(list(conf.GetAtomPosition(i)))
    s = []
    for i in ri:
        a = m.GetAtomWithIdx(i)
        me = [x.GetIdx() for x in a.GetNeighbors() if x.GetSymbol() == 'C' and x.GetIdx() not in ri]
        if not me:
            continue
        h = [x.GetIdx() for x in a.GetNeighbors() if x.GetSymbol() == 'H'][0]
        s.append(np.sign(np.dot(pos(me[0]) - pos(h), n)))   # メチル基が H より n 側にあるか
    return 'cis' if s[0] == s[1] else 'trans'

def chiral(smi):
    mir = smi.replace('@@', '!').replace('@', '@@').replace('!', '@')
    return Chem.CanonSmiles(smi) != Chem.CanonSmiles(mir)

cands = ['C[C@H]1CCC[C@@H](C)C1', 'C[C@H]1CCC[C@H](C)C1', 'C[C@H]1CC[C@@H](C)CC1', 'C[C@H]1CC[C@H](C)CC1']
res = {}
for s in cands:
    pos = '1,3' if 'CCC[' in s else '1,4'
    res[(pos, cis_trans(s))] = chiral(s)
print(res)
assert res == {('1,3', 'cis'): False, ('1,3', 'trans'): True, ('1,4', 'cis'): False, ('1,4', 'trans'): False}

# モデルのエネルギー：同じ側の関係 1,3 → (e,e)/(a,a)，1,4 → (a,e)/(e,a)
def confs(pos, ct):
    same_ae = {'1,3': False, '1,4': True}[pos]
    mixed = same_ae if ct == 'cis' else not same_ae
    return [('a', 'e'), ('e', 'a')] if mixed else [('e', 'e'), ('a', 'a')]
def E(pos, c):
    e = 8 * c.count('a')
    if pos == '1,3' and c == ('a', 'a'): e += 15
    return e
tab = {(p, ct): sorted(E(p, c) for c in confs(p, ct)) for p in ('1,3', '1,4') for ct in ('cis', 'trans')}
print(tab)
assert tab == {('1,3', 'cis'): [0, 31], ('1,3', 'trans'): [8, 8], ('1,4', 'cis'): [8, 8], ('1,4', 'trans'): [0, 16]}
r = 10 ** (-16000 / (2.30 * 8.31 * 300))
assert f'{r:.2g}' == '0.0016', r
print('ORG-STR-054 OK', r)
