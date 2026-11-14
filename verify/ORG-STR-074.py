# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-STR-074/item.json'
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
from itertools import product
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')

# 鎖式 C4H8O（重原子 C4 O1，二重結合 1 つ）をプリューファー列で総当たり
def prufer_edges(seq, n):
    deg = [1]*n
    for s in seq: deg[s] += 1
    edges = []
    for s in seq:
        leaf = min(i for i in range(n) if deg[i] == 1)
        edges.append((leaf, s)); deg[leaf] -= 1; deg[s] -= 1
    u, v = [i for i in range(n) if deg[i] == 1]
    edges.append((u, v)); return edges
iso = set()
for seq in product(range(5), repeat=3):
    E = prufer_edges(seq, 5)
    for o in range(5):
        for dbl in range(4):
            rw = Chem.RWMol()
            for i in range(5): rw.AddAtom(Chem.Atom(8 if i == o else 6))
            for k, (a, b) in enumerate(E):
                rw.AddBond(a, b, Chem.BondType.DOUBLE if k == dbl else Chem.BondType.SINGLE)
            try:
                m = rw.GetMol(); Chem.SanitizeMol(m)
            except Exception:
                continue
            if CalcMolFormula(m) == 'C4H8O':
                iso.add(Chem.MolToSmiles(m))
P = lambda s: Chem.MolFromSmarts(s)
ketones = sorted(s for s in iso if Chem.MolFromSmiles(s).HasSubstructMatch(P('[#6][CX3](=O)[#6]')))
aldeh = sorted(s for s in iso if Chem.MolFromSmiles(s).HasSubstructMatch(P('[CX3H1](=O)')))
alc = sorted(s for s in iso if Chem.MolFromSmiles(s).HasSubstructMatch(P('[CX4][OX2H]')))
enol = sorted(s for s in iso if Chem.MolFromSmiles(s).HasSubstructMatch(P('[CX3]=[CX3][OX2H]')))
print('鎖式C4H8O', len(iso)); print('ケトン', ketones); print('アルデヒド', aldeh); print('アルコール', alc); print('エノール(除外)', enol)
assert ketones == [Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)CC'))]
assert set(aldeh) == {Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in ('CCCC=O', 'CC(C)C=O')}
assert set(alc) == {Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in ('C=CCCO', 'CC=CCO', 'C=CC(C)O', 'C=C(C)CO')}
chiral = [s for s in alc if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True)]
print('不斉炭素をもつアルコール', chiral)
assert chiral == [Chem.MolToSmiles(Chem.MolFromSmiles('C=CC(C)O'))]
# 脱水（OH と隣の C 上の H）：2-メチル-1-プロパノール → 2-メチルプロペンだけ，1-ブタノール → 1-ブテン
def dehyd(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; c = o.GetNeighbors()[0]
    for n in c.GetNeighbors():
        if n.GetSymbol() == 'C' and n.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(o.GetIdx()); p = rw.GetMol()
            for a in p.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
            Chem.SanitizeMol(p); out.add(Chem.MolToSmiles(p))
    return out
assert dehyd('CC(C)CO') == {Chem.MolToSmiles(Chem.MolFromSmiles('C=C(C)C'))}
assert dehyd('CCCCO') == {Chem.MolToSmiles(Chem.MolFromSmiles('C=CCC'))}
# ヨードホルム：2-ブタノンは CH3CO- をもつ
assert Chem.MolFromSmiles('CC(=O)CC').HasSubstructMatch(P('[CH3]C(=O)'))
assert not Chem.MolFromSmiles('CC(C)C=O').HasSubstructMatch(P('[CH3]C(=O)'))
# エチルビニルエーテルと 3-ブテン-1-オールの OH の有無
assert not Chem.MolFromSmiles('C=COCC').HasSubstructMatch(P('[OX2H]')) and Chem.MolFromSmiles('C=CCCO').HasSubstructMatch(P('[OX2H]'))
for s in ('C=COCC', 'C=CCCO'): assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H8O'
print('OK')
