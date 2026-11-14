# 検証スクリプト ORG-ALC-015  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ALC-015/item.json'
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
alk = {1: ['C'], 2: ['CC'], 3: ['CCC', 'C(C)C'], 4: ['CCCC', 'C(C)CC', 'CC(C)C', 'C(C)(C)C']}
ethers = set()
for a in range(1, 5):
    b = 5 - a
    for x in alk[a]:
        for y in alk[b]:
            m = Chem.MolFromSmiles('O(' + x + ')' + y)
            ethers.add(Chem.MolToSmiles(m))
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C5H12O' for s in ethers)
assert len(ethers) == 6, ethers
chiral = [s for s in ethers if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True)]
assert chiral == [Chem.MolToSmiles(Chem.MolFromSmiles('COC(C)CC'))], chiral
# 2 種類のアルコール R1OH, R2OH の分子間脱水
alcs = ['CC', 'C(C)C']
prods = {Chem.MolToSmiles(Chem.MolFromSmiles('O(' + p + ')' + q)) for p in alcs for q in alcs}
assert len(prods) == 3
print('エーテル', sorted(ethers), '分子間脱水', prods)

print('ALL OK')
