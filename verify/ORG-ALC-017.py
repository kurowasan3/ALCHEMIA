# 検証スクリプト ORG-ALC-017  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-ALC-017/item.json'
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

import sympy as sp
from rdkit import Chem
x = sp.Symbol('x')
sols = sorted(sp.solve(sp.Eq(2*x*(1 - x), sp.Rational(42, 100)), x))
assert sols == [sp.Rational(3, 10), sp.Rational(7, 10)], sols
M_chi3 = 12 + 1.0 + 3*127
assert M_chi3 == 394
frac = round(35.9 / M_chi3, 4) / 0.100
pred = {float(s): float(1 - (1 - s)**2) for s in sols}
best = min(pred, key=lambda s: abs(pred[s] - frac))
assert best == 0.7 and abs(pred[0.7] - 0.91) < 1e-9 and abs(pred[0.3] - 0.51) < 1e-9
Ma = 4*12 + 6*1.0 + 4*16 + 40; Mp = 6*12 + 10*1.0 + 4*16 + 40
assert (Ma, Mp) == (158, 186)
pct = 7*Ma / (7*Ma + 3*Mp) * 100
assert sig(pct, 3) == 66.5, pct
f = 2*x*(1 - x)
xm = sp.solve(sp.diff(f, x), x)[0]
assert xm == sp.Rational(1, 2) and f.subs(x, xm) == sp.Rational(1, 2)
mk = Chem.MolFromSmarts('[CH3]C(=O)[#6]')
assert [Chem.MolFromSmiles(s).HasSubstructMatch(mk) for s in ('CC(C)=O', 'CCC(C)=O', 'CCC(=O)CC')] == [True, True, False]
print('x 候補', sols, 'メチルケトン割合', frac, '質量%', pct)

print('ALL OK')
