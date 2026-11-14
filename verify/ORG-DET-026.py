# 検証スクリプト ORG-DET-026  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-DET-026/item.json'
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
for e in ['CH3COCH3 + 3I2 + 4NaOH -> CH3COONa + CHI3 + 3NaI + 3H2O',
          'C3H8O + 4I2 + 6NaOH -> CH3COONa + CHI3 + 5NaI + 5H2O',
          'C6H10O2 + 6I2 + 8NaOH -> C4H4O4Na2 + 2CHI3 + 6NaI + 6H2O']:
    assert check_equation(e) is None, e
assert CalcMolFormula(Chem.MolFromSmiles('CC(=O)CCC(C)=O')) == 'C6H10O2'
CHI3 = 12 + 1.0 + 3*127; I2 = 254
assert CHI3 == 394
assert abs(3.94/CHI3 - 0.0100) < 1e-9
assert round(0.0100*I2, 2) == 2.54
nZ = 1.14/114; nI2 = 15.2/I2; nC = 7.88/CHI3
assert abs(nZ - 0.0100) < 1e-9 and abs(nC/nZ - 2) < 1e-9 and abs(nI2/nC - 3) < 0.02
print('ヨードホルム反応の量的関係 OK')

print('ALL OK')
