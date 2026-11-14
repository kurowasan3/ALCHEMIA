# 検証スクリプト ORG-HYD-029  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-HYD-029/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-HYD-029：問題文のモデルを RDKit で実装して全アルケンに適用 ----
from rdkit import Chem
def cdeg(m, i): return sum(1 for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == 'C')
def hcount(m, i): return 4 - sum(int(b.GetBondTypeAsDouble()) for b in m.GetAtomWithIdx(i).GetBonds())
def products(alkene, X=17):
    m0 = Chem.MolFromSmiles(alkene); cats = []
    for b in m0.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE: continue
        for hp, cat in ((b.GetBeginAtomIdx(), b.GetEndAtomIdx()), (b.GetEndAtomIdx(), b.GetBeginAtomIdx())):
            rw = Chem.RWMol(m0); rw.GetBondBetweenAtoms(hp, cat).SetBondType(Chem.BondType.SINGLE); Chem.SanitizeMol(rw)
            cats.append((cdeg(rw, cat), rw, cat))
    best = max(c[0] for c in cats); finals = []
    for d, rw, cat in cats:
        if d != best: continue
        cand = []
        for n in [x.GetIdx() for x in rw.GetAtomWithIdx(cat).GetNeighbors() if x.GetSymbol() == 'C']:
            if hcount(rw, n) >= 1 and cdeg(rw, n) > d: cand.append((cdeg(rw, n), rw, n, 'H'))
            for mth in [x.GetIdx() for x in rw.GetAtomWithIdx(n).GetNeighbors() if x.GetSymbol() == 'C' and x.GetIdx() != cat and x.GetDegree() == 1]:
                r2 = Chem.RWMol(rw); r2.RemoveBond(n, mth); r2.AddBond(cat, mth, Chem.BondType.SINGLE); Chem.SanitizeMol(r2)
                if cdeg(r2, n) > d: cand.append((cdeg(r2, n), r2, n, 'Me'))
        if cand:
            top = max(c[0] for c in cand); finals += [(c[1], c[2], c[3]) for c in cand if c[0] == top]
        else: finals.append((rw, cat, '-'))
    out = set()
    for rw, cat, how in finals:
        r = Chem.RWMol(rw); i = r.AddAtom(Chem.Atom(X)); r.AddBond(cat, i, Chem.BondType.SINGLE); Chem.SanitizeMol(r)
        out.add((Chem.MolToSmiles(r), how))
    return out
can = Chem.CanonSmiles
def nchiral(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=True))
# (1)(2)
p1 = products('C=CC(C)C'); p2 = products('C=CC(C)(C)C'); print(p1, p2)
assert p1 == {(can('CCC(C)(C)Cl'), 'H')} and p2 == {(can('CC(C)C(C)(C)Cl'), 'Me')}
assert nchiral('CCC(C)(C)Cl') == 0 and nchiral('CC(C)C(C)(C)Cl') == 0
assert nchiral('CC(C)C(C)Cl') == 1 and nchiral('CC(C)(C)C(C)Cl') == 1
# 鎖式アルケンの全列挙
def alkenes(skels):
    S = set()
    for s in skels:
        m = Chem.MolFromSmiles(s)
        for b in m.GetBonds():
            if b.GetBeginAtom().GetTotalNumHs() and b.GetEndAtom().GetTotalNumHs():
                rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE); Chem.SanitizeMol(rw); S.add(Chem.MolToSmiles(rw))
    return S
A5 = alkenes(['CCCCC', 'CC(C)CC', 'CC(C)(C)C']); print('C5H10', len(A5)); assert len(A5) == 5
t5 = can('CCC(C)(C)Cl')
hit5 = {a for a in A5 if {x for x, _ in products(a)} == {t5}}; print(hit5)
assert hit5 == {can('C=C(C)CC'), can('CC=C(C)C'), can('C=CC(C)C')}
assert {x for x,_ in products('CC=CCC')} == {can('CCCC(C)Cl'), can('CCC(Cl)CC')}
A6 = alkenes(['CCCCCC', 'CC(C)CCC', 'CCC(C)CC', 'CC(C)(C)CC', 'CC(C)C(C)C']); print('C6H12', len(A6)); assert len(A6) == 13
t6 = can('CC(C)C(C)(C)Cl')
hit6 = {a for a in A6 if {x for x, _ in products(a)} == {t6}}; print(hit6)
assert hit6 == {can('C=C(C)C(C)C'), can('CC(C)=C(C)C'), can('C=CC(C)(C)C')}
assert products('C=CC(C)CC') == {(can('CCC(C)(Cl)CC'), 'H')}
# (5) 水の付加
p5 = products('C=CC(C)(C)C', X=8); print(p5); assert {x for x,_ in p5} == {can('CC(C)C(C)(C)O')}
m = Chem.MolFromSmiles('CC(C)C(C)(C)O')
cOH = [a for a in m.GetAtoms() if a.GetSymbol()=='O'][0].GetNeighbors()[0]
assert cOH.GetTotalNumHs() == 0   # 第三級アルコール
print('ALL OK')
