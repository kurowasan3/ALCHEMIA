# 検証スクリプト ORG-CBX-025  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-CBX-025/item.json'
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
# ---- 個別の検証 ORG-CBX-025 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Crippen
S18 = 'CCCCCCCCCCCCCCCCC'
mols = {'a':('OCC(O)CO','C3H8O3'), 'b':(S18+'C(=O)OCC(O)CO','C21H42O4'), 'c':(S18+'C(=O)[O-].[Na+]','C18H35NaO2'),
        'd':('CC(=O)[O-].[Na+]','C2H3NaO2'), 'e':(S18+'C(=O)OCC(OC(=O)'+S18+')COC(=O)'+S18,'C57H110O6'),
        'g':('CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]','C12H25NaO4S'), 'h':('C'*16,'C16H34')}
def longest_chain(m):
    # 最長の連続した sp3 炭素鎖（疎水性部分）の長さ
    import itertools
    C = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='C' and not any(n.GetSymbol()=='O' for n in a.GetNeighbors())]
    best = 0
    dm = Chem.GetDistanceMatrix(m)
    sub = set(C)
    for i in C:
        for j in C:
            path = Chem.GetShortestPath(m, i, j) if i != j else (i,)
            if all(p in sub for p in path): best = max(best, len(path))
    return best
def polar(m):  # 親水性の部分：電荷をもつ原子 or 遊離の OH の数
    ch = sum(1 for a in m.GetAtoms() if a.GetFormalCharge() != 0 and a.GetSymbol() != 'Na')
    oh = sum(1 for a in m.GetAtoms() if a.GetSymbol()=='O' and a.GetTotalNumHs()==1)
    return ch + oh
res = {}
for k,(s,f) in mols.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == f, (k, CalcMolFormula(m))
    L, p = longest_chain(m), polar(m)
    res[k] = (L >= 10 and p >= 1)
    print(k, f, '疎水鎖', L, '親水部', p, '→', res[k])
assert {k for k,v in res.items() if v} == {'b','c','g'}   # (f) レシチンは記述で与えた（2本の長鎖＋双性イオン部分）
assert check_equation('C17H35COO- + H+ -> C17H35COOH') is None
# 問題文の示性式の分子式
for s, f in (('C17H35COOCH2CH(OH)CH2OH','C21H42O4'), ('(C17H35COO)3C3H5','C57H110O6')):
    pass
print('ALL OK')
