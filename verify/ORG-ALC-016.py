# 検証スクリプト ORG-ALC-016  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ALC-016/item.json'
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
prods = {'a': 'CCCO', 'b': 'CC(O)CC', 'c': 'CC(C)(O)CC', '3-pentanol': 'CCC(O)CC'}
def deg(smi):
    m = Chem.MolFromSmiles(smi)
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c = o.GetNeighbors()[0]
    return sum(1 for x in c.GetNeighbors() if x.GetSymbol() == 'C')
assert [deg(prods[k]) for k in 'abc'] == [1, 2, 3]
assert [CalcMolFormula(Chem.MolFromSmiles(prods[k])) for k in 'abc'] == ['C3H8O', 'C4H10O', 'C5H12O']
iod = Chem.MolFromSmarts('[CH3][CX4H1][OX2H1]')
assert [Chem.MolFromSmiles(prods[k]).HasSubstructMatch(iod) for k in 'abc'] == [False, True, False]
ch = lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
assert [ch(prods[k]) for k in 'abc'] == [0, 1, 0]
# 逆合成：プロパナール + C2H5 → 3-ペンタノール（カルボニル炭素に Et を付ける）
assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(O)CC')) == Chem.MolToSmiles(Chem.MolFromSmiles('CCC(CC)O'))
assert check_equation('C2H5MgBr + H2O -> C2H6 + MgBr(OH)') is None
print('生成物の級数・ヨードホルム・不斉 OK')

print('ALL OK')
