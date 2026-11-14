# 検証スクリプト ORG-STR-025  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-025/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b) and (' + ' in b or re.search(r'(^|\s)\d+\s*[A-Z(\[]', b)):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(nums) == n_items or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
def acids(core, k):
    m0 = Chem.AddHs(Chem.MolFromSmiles(core))
    from itertools import combinations
    hs = [a.GetIdx() for a in m0.GetAtoms() if a.GetSymbol() == 'H']
    res = set()
    for comb in combinations(hs, k):
        rw = Chem.RWMol(m0)
        for h in sorted(comb, reverse=True):
            nb = rw.GetAtomWithIdx(h).GetNeighbors()[0].GetIdx()
            rw.RemoveAtom(h)
            c = rw.AddAtom(Chem.Atom(6)); o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(nb, c, Chem.BondType.SINGLE); rw.AddBond(c, o1, Chem.BondType.DOUBLE); rw.AddBond(c, o2, Chem.BondType.SINGLE)
        mol = Chem.RemoveHs(rw.GetMol())
        Chem.SanitizeMol(mol)
        for iso in EnumerateStereoisomers(mol, options=StereoEnumerationOptions(onlyUnassigned=True, unique=True)):
            res.add(Chem.MolToSmiles(iso))
    return res
cc = acids('C=CC', 1)
ring = acids('C1CC1', 1)
assert len(cc) == 4, cc
assert len(ring) == 1, ring
for s in cc | ring:
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H6O2'
di = acids('CC', 2)
assert len(di) == 2 and all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H6O4' for s in di), di
unsat = lambda n, m: (2*n + 2 - m) // 2
assert unsat(6, 10) == 2 and unsat(8, 14) == 2 and unsat(4, 6) == 2
# (5) 水素付加で直鎖 & シス-トランスなし
def hydro(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C':
            b.SetBondType(Chem.BondType.SINGLE); b.SetStereo(Chem.BondStereo.STEREONONE)
    mm = m.GetMol(); Chem.SanitizeMol(mm)
    return Chem.MolToSmiles(mm, isomericSmiles=False)
butanoic = Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(=O)O'))
cand = {Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False) for s in cc if hydro(s) == butanoic}
assert cand == {Chem.MolToSmiles(Chem.MolFromSmiles('C=CCC(=O)O')), Chem.MolToSmiles(Chem.MolFromSmiles('CC=CC(=O)O'))}, cand
no_ct = [s for s in cand if '/' not in Chem.MolToSmiles(Chem.MolFromSmiles(s)) and s != Chem.MolToSmiles(Chem.MolFromSmiles('CC=CC(=O)O'))]
assert no_ct == [Chem.MolToSmiles(Chem.MolFromSmiles('C=CCC(=O)O'))]
print('C=C の酸', sorted(cc), '環の酸', ring, '2価の酸', di)

print('ALL OK')
