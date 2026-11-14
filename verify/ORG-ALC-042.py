# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-ALC-042/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from itertools import product
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
MASS = {'C': 12, 'H': 1, 'O': 16}
def nominal(m):
    return sum(MASS[a.GetSymbol()] + a.GetTotalNumHs()*1 for a in m.GetAtoms())
# 鎖式 C6H12O のケトンを総当たり（重原子7の木＋C=O）
def prufer_edges(seq, n):
    deg = [1]*n
    for s in seq: deg[s] += 1
    E = []
    for s in seq:
        leaf = min(i for i in range(n) if deg[i] == 1); E.append((leaf, s)); deg[leaf] -= 1; deg[s] -= 1
    u, v = [i for i in range(n) if deg[i] == 1]; E.append((u, v)); return E
kets = set()
for seq in product(range(7), repeat=5):
    E = prufer_edges(seq, 7)
    for o in range(7):
        oe = [k for k, (a, b) in enumerate(E) if o in (a, b)]
        if len(oe) != 1: continue
        rw = Chem.RWMol()
        for i in range(7): rw.AddAtom(Chem.Atom(8 if i == o else 6))
        for k, (a, b) in enumerate(E): rw.AddBond(a, b, Chem.BondType.DOUBLE if k == oe[0] else Chem.BondType.SINGLE)
        try:
            m = rw.GetMol(); Chem.SanitizeMol(m)
        except Exception: continue
        if CalcMolFormula(m) == 'C6H12O' and m.HasSubstructMatch(Chem.MolFromSmarts('[#6][CX3](=O)[#6]')):
            kets.add(Chem.MolToSmiles(m))
print('ケトン', len(kets), sorted(kets)); assert len(kets) == 6
def acyl(smi):
    m = Chem.MolFromSmiles(smi)
    c = m.GetSubstructMatch(Chem.MolFromSmarts('[#6][CX3](=O)[#6]'))
    res = []
    for side in (c[0], c[3]):
        other = c[3] if side == c[0] else c[0]
        fm = Chem.FragmentOnBonds(m, [m.GetBondBetweenAtoms(c[1], other).GetIdx()], addDummies=False)
        for p in Chem.GetMolFrags(fm, asMols=True):
            if p.HasSubstructMatch(Chem.MolFromSmarts('C=O')):
                res.append(nominal(p) - 1)   # 切断で付いた H を除いたアシリウムイオン
    return sorted(res)
tab = {s: acyl(s) for s in kets}
for s, v in tab.items(): print(s, v)
K_cand = sorted(s for s, v in tab.items() if v == [57, 71])
L_cand = sorted(s for s, v in tab.items() if v == [43, 85])
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert set(K_cand) == {can('CCC(=O)CCC'), can('CCC(=O)C(C)C')} and len(L_cand) == 4
ch3co = Chem.MolFromSmarts('[CH3]C(=O)')
assert not any(Chem.MolFromSmiles(s).HasSubstructMatch(ch3co) for s in K_cand)
assert all(Chem.MolFromSmiles(s).HasSubstructMatch(ch3co) for s in L_cand)
def alkenes(ket):
    alc = Chem.RWMol(Chem.MolFromSmiles(ket))
    o = [a for a in alc.GetAtoms() if a.GetSymbol() == 'O'][0]
    alc.GetBondBetweenAtoms(o.GetIdx(), o.GetNeighbors()[0].GetIdx()).SetBondType(Chem.BondType.SINGLE)
    m = alc.GetMol(); Chem.SanitizeMol(m)
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; c = o.GetNeighbors()[0]
    out = set()
    for n in c.GetNeighbors():
        if n.GetSymbol() == 'C' and n.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(o.GetIdx()); p = rw.GetMol()
            for a in p.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
            Chem.SanitizeMol(p)
            for st in EnumerateStereoisomers(p): out.add(Chem.MolToSmiles(st))
    return out
n = {s: len(alkenes(s)) for s in K_cand}; print(n)
assert n[can('CCC(=O)CCC')] == 4 and n[can('CCC(=O)C(C)C')] == 3
print('OK')
