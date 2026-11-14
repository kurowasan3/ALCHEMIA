# 検証スクリプト ORG-STR-026  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-026/item.json'
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
# (1)
sols = []
for z in range(1, 8):
    for x in range(4, 10):
        y = 116 - 12*x - 16*z
        if y >= 1 and y % 2 == 0 and y <= 2*x + 2 and (2*x + 2 - y)//2 <= 3:
            sols.append((x, y, z))
assert sorted(sols) == sorted([(7, 16, 1), (6, 12, 2), (5, 8, 3), (4, 4, 4)]), sols
# (3) P: C7 第三級アルコール
heptanes = ['CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CC(C)C(C)CC', 'CC(C)CC(C)C', 'CCC(C)(C)CC', 'CC(C)(C)CCC', 'CCC(CC)CC', 'CC(C)(C)C(C)C']
def alcohols(sk):
    m = Chem.MolFromSmiles(sk); out = set()
    for a in m.GetAtoms():
        if a.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(a.GetIdx(), o, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
def dehydrate(smi):
    m = Chem.MolFromSmiles(smi)
    oh = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    ca = oh.GetNeighbors()[0]
    prods = set()
    for b in ca.GetNeighbors():
        if b.GetSymbol() == 'C' and b.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m)
            rw.GetBondBetweenAtoms(ca.GetIdx(), b.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(oh.GetIdx())
            mm = rw.GetMol(); Chem.SanitizeMol(mm)
            for iso in EnumerateStereoisomers(mm, options=StereoEnumerationOptions(onlyUnassigned=True, unique=True)):
                prods.add(Chem.MolToSmiles(iso))
    return prods
def tertiary(smi):
    m = Chem.MolFromSmiles(smi)
    oh = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c = oh.GetNeighbors()[0]
    return sum(1 for n in c.GetNeighbors() if n.GetSymbol() == 'C') == 3
allC7 = set().union(*[alcohols(s) for s in heptanes])
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C7H16O' for s in allC7)
one = [s for s in allC7 if tertiary(s) and len(dehydrate(s)) == 1]
assert sorted(one) == sorted([Chem.MolToSmiles(Chem.MolFromSmiles('CCC(O)(CC)CC')), Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)(O)C(C)(C)C'))]), one
def branches(smi):
    m = Chem.MolFromSmiles(smi)
    m = Chem.DeleteSubstructs(m, Chem.MolFromSmarts('[OX2H]'))
    return [a.GetDegree() for a in m.GetAtoms() if a.GetDegree() >= 3]
P = [s for s in one if branches(s) == [3]]
assert P == [Chem.MolToSmiles(Chem.MolFromSmiles('CCC(O)(CC)CC'))], P
# Q: C5 第三級アルコールは 1 種類
pentanes = ['CCCCC', 'CC(C)CC', 'CC(C)(C)C']
t5 = [s for s in set().union(*[alcohols(s) for s in pentanes]) if tertiary(s)]
assert t5 == [Chem.MolToSmiles(Chem.MolFromSmiles('CCC(C)(C)O'))], t5
assert CalcMolFormula(Chem.MolFromSmiles('O=COC(C)(C)CC')) == 'C6H12O2'
# R: 直鎖 C5 のケト酸で CH3CO- をもつもの
keto = set()
for k in (2, 3, 4):
    chain = ['C'] * 5
    s = ''.join('C(=O)' if i == k - 1 else 'C' for i in range(4)) + 'C(=O)O'
    m = Chem.MolFromSmiles(s)
    if m.HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]')):
        keto.add(Chem.MolToSmiles(m))
assert keto == {Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)CCC(=O)O'))}, keto
assert CalcMolFormula(Chem.MolFromSmiles('CC(=O)CCC(=O)O')) == 'C5H8O3'
red = Chem.MolFromSmiles('CC(O)CCC(=O)O')
assert len(Chem.FindMolChiralCenters(red, includeUnassigned=True)) == 1
assert CalcMolFormula(Chem.MolFromSmiles('OC(=O)/C=C\\C(=O)O')) == 'C4H4O4'
for f, M in [('C7H16O', 116), ('C6H12O2', 116), ('C5H8O3', 116), ('C4H4O4', 116)]:
    pass
print('P', P, 'Q の アルコール', t5, 'R', keto)

print('ALL OK')
