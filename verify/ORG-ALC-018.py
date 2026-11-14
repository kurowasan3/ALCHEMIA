# 検証スクリプト ORG-ALC-018  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-ALC-018/item.json'
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
# 候補エステル（アルコール部分が CH3CH(OH)- をもつ第二級，酸化でケトン）
cands = {'iPr propanoate': ('CCC(=O)OC(C)C', 'CC', 'CC(C)O'),
         'sBu acetate': ('CC(=O)OC(C)CC', 'C', 'CC(O)CC'),
         '2-pentyl formate': ('O=COC(C)CCC', '', 'CC(O)CCC'),
         '3-me-2-butyl formate': ('O=COC(C)C(C)C', '', 'CC(O)C(C)C')}
for k, (e, r1, g) in cands.items():
    assert CalcMolFormula(Chem.MolFromSmiles(e)) == 'C6H12O2', k
def grignard_b(r1):   # R'C(OH)(CH3)2
    return can(f'{r1}C(C)(C)O') if r1 else can('CC(C)O')
def deg(smi):
    m = Chem.MolFromSmiles(smi); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    return sum(1 for x in o.GetNeighbors()[0].GetNeighbors() if x.GetSymbol() == 'C')
def alkenes(smi):
    m = Chem.MolFromSmiles(smi); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; c = o.GetNeighbors()[0]
    out = set()
    for nb in c.GetNeighbors():
        if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.RemoveBond(c.GetIdx(), o.GetIdx()); rw.RemoveAtom(o.GetIdx())
            rw.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            x = rw.GetMol(); Chem.SanitizeMol(x)
            for st in EnumerateStereoisomers(x): out.add(Chem.MolToSmiles(st))
    return out
res = {k: (deg(grignard_b(r1)), len(alkenes(grignard_b(r1)))) for k, (e, r1, g) in cands.items()}
ok = [k for k, (d, n) in res.items() if d == 3 and n == 2]
assert ok == ['iPr propanoate'], res
assert alkenes('CCC(C)(C)O') == {can('C=C(C)CC'), can('CC=C(C)C')}
assert grignard_b('') == can('CC(C)O') and deg(grignard_b('')) == 2   # ギ酸エステル → 第二級
M = 6*12 + 12*1.0 + 2*16; assert M == 116
assert round((0.500 - 2*11.6/M)*22.4, 2) == 6.72
for e in ['CH3MgBr + H2O -> CH4 + MgBr(OH)']: assert check_equation(e) is None
# (4) (a) の組：ケトン + RMgBr で F
rxn = __import__('rdkit.Chem.AllChem', fromlist=['x']).ReactionFromSmarts('[C:1]=[O:2].[CH3,CH2:3][Mg]>>[C:1]([*:3])[O:2]')
def addR(carbonyl, rmg):
    p = rxn.RunReactants((Chem.MolFromSmiles(carbonyl), Chem.MolFromSmiles(rmg)))[0][0]; Chem.SanitizeMol(p); return Chem.MolToSmiles(p)
assert addR('CC(C)=O', 'CC[Mg]') == can('CC(C)(O)CC') and addR('CC(=O)CC', 'C[Mg]') == can('CC(C)(O)CC')
print('E・F・G の決定，CH4 6.72 L OK')

print('ALL OK')
