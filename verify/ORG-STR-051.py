# 検証スクリプト ORG-STR-051（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問の数', np_, na)
    return it

# ---- 異性体生成器（単独で動くよう埋め込み） ----
# 小さな構造異性体生成器: C と O の数を指定して、全ての連結な分子（H は暗黙）を列挙する
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
VAL = {'C': 4, 'O': 2}
def enumerate_isomers(nC, nO, maxbond=3):
    start = set()
    for a in ('C', 'O'):
        start.add(Chem.CanonSmiles(a))
    seen = set(start); frontier = set(start)
    while frontier:
        new = set()
        for smi in frontier:
            m = Chem.MolFromSmiles(smi)
            cC = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'C'); cO = m.GetNumAtoms() - cC
            rw0 = Chem.RWMol(m); Chem.Kekulize(rw0, clearAromaticFlags=True)
            cands = []
            # 原子の追加
            for el, ok in (('C', cC < nC), ('O', cO < nO)):
                if not ok: continue
                for a in rw0.GetAtoms():
                    for bo in (1, 2, 3):
                        rw = Chem.RWMol(rw0); i = rw.AddAtom(Chem.Atom(el))
                        rw.AddBond(a.GetIdx(), i, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[bo])
                        cands.append(rw)
            # 結合の追加・次数の増加
            n = rw0.GetNumAtoms()
            for i in range(n):
                for j in range(i + 1, n):
                    b = rw0.GetBondBetweenAtoms(i, j)
                    rw = Chem.RWMol(rw0)
                    if b is None:
                        rw.AddBond(i, j, Chem.BondType.SINGLE)
                    else:
                        o = b.GetBondTypeAsDouble()
                        if o >= maxbond: continue
                        rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE if o == 1 else Chem.BondType.TRIPLE)
                    cands.append(rw)
            for rw in cands:
                ok = True
                for a in rw.GetAtoms():
                    if sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) > VAL[a.GetSymbol()]: ok = False; break
                if not ok: continue
                try:
                    mm = rw.GetMol(); 
                    for a in mm.GetAtoms(): a.SetNoImplicit(False)
                    Chem.SanitizeMol(mm)
                    s = Chem.MolToSmiles(mm)
                except Exception:
                    continue
                if s not in seen:
                    seen.add(s); new.add(s)
        frontier = new
    out = []
    for s in seen:
        m = Chem.MolFromSmiles(s)
        cC = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'C'); cO = m.GetNumAtoms() - cC
        if cC == nC and cO == nO: out.append(s)
    return out

it = check_item('/home/claude/chem6/outh/items/ORG-STR-051/item.json')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
assert len([s for s in enumerate_isomers(4, 1) if F(Chem.MolFromSmiles(s)) == 'C4H10O']) == 7
L = [s for s in enumerate_isomers(5, 2) if F(Chem.MolFromSmiles(s)) == 'C5H8O2']
print('C5H8O2 の構造', len(L))
mol = Chem.MolFromSmiles
ester = Chem.MolFromSmarts('[#6,#1][CX3](=O)[OX2][#6]')
OO = Chem.MolFromSmarts('[#8]~[#8]'); unsat = Chem.MolFromSmarts('[#6]=,#[#6]')
def is_ester(m):
    # C(=O)-O-C で、O の相手の C がカルボニル C の隣（三員環）でもよいものを拾う
    patt = Chem.MolFromSmarts('[CX3](=O)[OX2][#6]')
    return m.HasSubstructMatch(patt)
E = [s for s in L if is_ester(mol(s)) and not mol(s).HasSubstructMatch(OO) and not mol(s).HasSubstructMatch(unsat)]
def ester_ring_size(m):
    patt = Chem.MolFromSmarts('[CX3](=O)[OX2][#6]')
    ri = m.GetRingInfo()
    for c, o1, o2, x in m.GetSubstructMatches(patt):
        b = m.GetBondBetweenAtoms(c, o2).GetIdx()
        if ri.NumBondRings(b):
            return min(len(r) for r in ri.BondRings() if b in r)
    return 0
three = [s for s in E if ester_ring_size(mol(s)) == 3]
E2 = [s for s in E if ester_ring_size(mol(s)) != 3]
carbo = [s for s in E2 if ester_ring_size(mol(s)) == 0]
lact = [s for s in E2 if ester_ring_size(mol(s)) > 0]
print('三員環ラクトン（除外）', len(three), three)
print('炭素環のエステル', len(carbo), sorted(carbo))
print('ラクトン', len(lact), sorted(lact), [ester_ring_size(mol(s)) for s in sorted(lact)])
assert len(carbo) == 6 and len(lact) == 9
assert all(mol(s).GetRingInfo().NumRings() == 1 for s in E2)
def nstereo(s): return len(set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(mol(s))))
st = {s: nstereo(s) for s in carbo}
print('炭素環エステルの立体異性体数', st)
assert [s for s in carbo if st[s] > 1] == [Chem.CanonSmiles('CC1CC1OC=O')] and st[Chem.CanonSmiles('CC1CC1OC=O')] == 4
# ラクトンを加水分解した Y（ヒドロキシ酸）を作る
from rdkit.Chem import AllChem
def hydrolyze(s):
    m = mol(s); patt = Chem.MolFromSmarts('[CX3](=O)[OX2][#6]')
    c, o1, o2, x = m.GetSubstructMatch(patt)
    rw = Chem.RWMol(m); rw.RemoveBond(c, o2)
    i = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, i, Chem.BondType.SINGLE)
    p = rw.GetMol(); 
    for a in p.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(p)
    return {Chem.MolToSmiles(p)}
Ys = {}
for s in lact:
    pr = hydrolyze(s); assert len(pr) == 1, (s, pr)
    y = pr.pop(); assert F(mol(y)) == 'C5H10O3'; Ys[y] = s
for s in carbo:
    pr = hydrolyze(s); 
    for p in pr: assert '.' in p, (s, p)   # 2 分子に分かれる
print('Y の候補', len(Ys))
def nchiral(s): return len(Chem.FindMolChiralCenters(mol(s), includeUnassigned=True, useLegacyImplementation=False))
ox1 = AllChem.ReactionFromSmarts('[CH2:1][OH:2]>>[CH1:1]=[O:2]')
ox2 = AllChem.ReactionFromSmarts('[CH1:1]([#6:3])([#6:4])[OH:2]>>[C:1]([#6:3])([#6:4])=[O:2]')
iodo = Chem.MolFromSmarts('[CH3][CX3](=O)[#6,#1]')
ok = []
for y in Ys:
    if nchiral(y) != 1: continue
    m = mol(y); Z = None
    for rx in (ox1, ox2):
        ps = rx.RunReactants((m,))
        if ps:
            z = ps[0][0]; Chem.SanitizeMol(z); Z = Chem.MolToSmiles(z)
    assert Z and F(mol(Z)) == 'C5H8O3'
    print(y, '→', Z, '不斉', nchiral(Z), 'CH3CO', mol(Z).HasSubstructMatch(iodo))
    if nchiral(Z) == 0 and mol(Z).HasSubstructMatch(iodo): ok.append(y)
print('条件を満たす Y', ok)
assert ok == [Chem.CanonSmiles('CC(O)CCC(=O)O')]
X = Ys[ok[0]]; print('X =', X, '立体異性体', nstereo(X)); assert nstereo(X) == 2
print('OK')
