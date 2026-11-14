# 検証スクリプト THE-RDX-026  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-RDX-026/item.json'
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
assert CalcMolFormula(Chem.MolFromSmiles('O=C1/C(=C2\\Nc3ccccc3C2=O)Nc2ccccc12')) == 'C16H10N2O2'
assert CalcMolFormula(Chem.MolFromSmiles('[Na]Oc1c(-c2[nH]c3ccccc3c2O[Na])[nH]c2ccccc12')) == 'C16H10N2Na2O2'
M_ind = 16*12 + 10*1.0 + 2*14 + 2*16
M_dt = 2*23 + 2*32 + 4*16
assert M_ind == 262 and M_dt == 174
n = 5.24 / M_ind
assert abs(n - 0.0200) < 1e-12
assert sig(n * M_dt, 3) == 3.48
assert sig(n / 2 * 22.4, 3) == 0.224
x = (-2 + 8) / 2
assert x == 3 and (-2 + 6) == 4
assert check_equation('2C16H10N2O2Na2 + O2 + 2H2O -> 2C16H10N2O2 + 4NaOH') is None
print('量的関係 OK')

print('ALL OK')
