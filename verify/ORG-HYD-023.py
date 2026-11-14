ITEM = '/home/claude/chem2/outh/items/ORG-HYD-023/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
import numpy as np
from rdkit.Geometry import Point3D
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
mirror = lambda s: can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
MESO = 'C[C@@H](Br)[C@@H](Br)C'; RAC = 'C[C@@H](Br)[C@H](Br)C'
assert mirror(MESO) == can(MESO) and mirror(RAC) != can(RAC)
def anti_add(alkene_smi, seed=1):
    m = Chem.AddHs(Chem.MolFromSmiles(alkene_smi)); assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(m); conf = m.GetConformer()
    db = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    i, j = db.GetBeginAtomIdx(), db.GetEndAtomIdx()
    P = lambda k: np.array(conf.GetAtomPosition(k))
    nb = [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() != j]
    n = np.cross(P(nb[0]) - P(i), P(j) - P(i)); n /= np.linalg.norm(n)
    rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.SINGLE)
    bi = rw.AddAtom(Chem.Atom(35)); bj = rw.AddAtom(Chem.Atom(35))
    rw.AddBond(i, bi, Chem.BondType.SINGLE); rw.AddBond(j, bj, Chem.BondType.SINGLE)
    c = rw.GetConformer()
    for k, pos in [(bi, P(i) + 1.9 * n), (bj, P(j) - 1.9 * n)]: c.SetAtomPosition(k, Point3D(*pos))
    mm = rw.GetMol(); Chem.SanitizeMol(mm); Chem.AssignStereochemistryFrom3D(mm)
    return can(Chem.MolToSmiles(Chem.RemoveHs(mm)))
for seed in (1, 5):
    assert anti_add('C/C=C/C', seed) == can(MESO)
    assert anti_add('C/C=C\\C', seed) in (can(RAC), mirror(RAC))
def geom(smi, l2, l3, angle, seed):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(m); conf = m.GetConformer()
    pick = lambda c, sym: [n.GetIdx() for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetSymbol() == sym][0]
    T.SetDihedralDeg(conf, pick(1, l2), 1, 3, pick(3, l3), angle)
    d = T.GetDihedralDeg(conf, 0, 1, 3, 5)
    return 'same' if abs(d) < 90 else 'opp'
for seed in (1, 7, 42):
    assert geom(MESO, 'Br', 'Br', 180, seed) == 'opp'    # (2) P -> トランス
    assert geom(RAC, 'Br', 'Br', 180, seed) == 'same'    # (2) Q -> シス
    assert geom(MESO, 'H', 'Br', 180, seed) == 'same' and geom(MESO, 'Br', 'H', 180, seed) == 'same'   # (3) P
    assert geom(RAC, 'H', 'Br', 180, seed) == 'opp' and geom(RAC, 'Br', 'H', 180, seed) == 'opp'       # (3) Q
    assert geom(MESO, 'Br', 'Br', 0, seed) == 'same'     # (4) 同じ向きなら P -> シス
print('立体の判定OK')
# 解答図のニューマン投影の角度（手前 C2：Br0, CH3 120, H 240）を検査
F = {'Br': 0, 'CH3': 120, 'H': 240}
backs = {'a': {'Br': 180, 'CH3': 300, 'H': 60}, 'b': {'Br': 180, 'CH3': 60, 'H': 300}, 'c': {'Br': 60, 'CH3': 180, 'H': 300}, 'd': {'Br': 60, 'CH3': 300, 'H': 180}}
inv = lambda b: all((F[k] + 180) % 360 == b[k] for k in F)          # 中心対称（メソ）
c2 = lambda b: all((180 - F[k]) % 360 == b[k] for k in F)           # C2 対称（キラル）
rot = lambda b, d: {k: (v + d) % 360 for k, v in b.items()}
assert inv(backs['a']) and c2(backs['b'])
assert backs['c'] == rot(backs['a'], -120) and backs['d'] == rot(backs['b'], -120)
ang = lambda x, y: min(abs(x - y) % 360, 360 - abs(x - y) % 360)
assert ang(F['CH3'], backs['a']['CH3']) == 180 and ang(F['CH3'], backs['b']['CH3']) == 60
assert ang(F['H'], backs['c']['Br']) == 180 and ang(F['CH3'], backs['c']['CH3']) == 60
assert ang(F['H'], backs['d']['Br']) == 180 and ang(F['CH3'], backs['d']['CH3']) == 180
svg = open(ITEM.replace('item.json', 'solution.svg'), encoding='utf-8').read(); assert 'viewBox' in svg and '<script' not in svg

print('全検証OK:', ITEM)
