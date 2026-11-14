# 検証スクリプト ORG-ARO-064（グループ g09）  python3 -I この.py
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

# ---- 構造の列挙（単独で動くよう埋め込み） ----
# 出発分子に原子を 1 個ずつ（鎖状に）付け加えて構造を列挙する（新しい環はつくらない）
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
VAL2 = {'C': 4, 'O': 2, 'N': 3}
BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
def grow(start, add, maxbond=2):
    """start: SMILES, add: dict 元素→個数"""
    from collections import Counter
    def cnt(m):
        return Counter(a.GetSymbol() for a in m.GetAtoms())
    base = cnt(Chem.MolFromSmiles(start))
    target = base + Counter(add)
    seen = {Chem.CanonSmiles(start)}; frontier = set(seen)
    while frontier:
        new = set()
        for smi in frontier:
            m = Chem.MolFromSmiles(smi); c = cnt(m)
            rw0 = Chem.RWMol(m); Chem.Kekulize(rw0, clearAromaticFlags=True)
            for el in add:
                if c[el] >= target[el]: continue
                for a in rw0.GetAtoms():
                    for bo in range(1, maxbond + 1):
                        rw = Chem.RWMol(rw0); i = rw.AddAtom(Chem.Atom(el)); rw.AddBond(a.GetIdx(), i, BT[bo])
                        if any(sum(b.GetBondTypeAsDouble() for b in x.GetBonds()) > VAL2[x.GetSymbol()] for x in rw.GetAtoms()): continue
                        try:
                            mm = rw.GetMol()
                            for x in mm.GetAtoms(): x.SetNoImplicit(False)
                            Chem.SanitizeMol(mm); s = Chem.MolToSmiles(mm)
                        except Exception:
                            continue
                        if s not in seen: seen.add(s); new.add(s)
            # 既存の結合の次数を上げる（環の外の結合のみ）
            for b in rw0.GetBonds():
                if b.IsInRing() or b.GetBondTypeAsDouble() >= maxbond: continue
                rw = Chem.RWMol(rw0); rw.GetBondWithIdx(b.GetIdx()).SetBondType(BT[int(b.GetBondTypeAsDouble()) + 1])
                if any(sum(bb.GetBondTypeAsDouble() for bb in x.GetBonds()) > VAL2[x.GetSymbol()] for x in rw.GetAtoms()): continue
                try:
                    mm = rw.GetMol()
                    for x in mm.GetAtoms(): x.SetNoImplicit(False)
                    Chem.SanitizeMol(mm); s = Chem.MolToSmiles(mm)
                except Exception:
                    continue
                if s not in seen: seen.add(s); new.add(s)
        frontier = new
    return [s for s in seen if cnt(Chem.MolFromSmiles(s)) == target]

it = check_item('/home/claude/chem6/out/items/ORG-ARO-064/item.json')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
mol = Chem.MolFromSmiles
L = [s for s in grow('c1ccccc1', {'C': 2, 'N': 1, 'O': 2}) if F(mol(s)) == 'C8H9NO2']
print('C8H9NO2（ベンゼン環＋鎖状部分）', len(L))
nh2 = Chem.MolFromSmarts('[NX3H2][#6;!$(C=O)]'); cooh = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
A = [s for s in L if len(mol(s).GetSubstructMatches(nh2)) == 1 and len(mol(s).GetSubstructMatches(cooh)) == 1]
print('(1)', len(A)); assert len(A) == 17
arN = Chem.MolFromSmarts('c[NH2]')
B = [s for s in A if mol(s).HasSubstructMatch(arN)]
print('(2)', len(B)); assert len(B) == 13
C = [s for s in A if Chem.FindMolChiralCenters(mol(s), includeUnassigned=True, useLegacyImplementation=False)]
print('(3)', C); assert C == [Chem.CanonSmiles('NC(C(=O)O)c1ccccc1')]
def ring_H_types(s):
    m = Chem.AddHs(mol(s)); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[h.GetIdx()] for a in m.GetAtoms() if a.GetIsAromatic() for h in a.GetNeighbors() if h.GetSymbol() == 'H'})
W = [s for s in A if s not in B and ring_H_types(s) == 2]
print('(4)', W, {s: ring_H_types(s) for s in A if s not in B}); assert W == [Chem.CanonSmiles('NCc1ccc(C(=O)O)cc1')]
print('OK')
