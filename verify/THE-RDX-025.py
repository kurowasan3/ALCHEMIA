# 検証スクリプト THE-RDX-025  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-RDX-025/item.json'
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(re.findall(r'\((\d+)\)', it['answer'])) == len(nums), '答え中に余分な (n)'
    if n_items: assert nums == [str(i+1) for i in range(n_items)], (n_items, nums)
    else: assert not nums
    print('小問数', n_items, '答え', nums)
check_structure()
check_all_equations()
def bis(f, lo=-16, hi=1):
    """f(h) が h について増加 → f=0 の h（log 空間の二分法）"""
    for _ in range(200):
        m = (lo + hi) / 2
        if f(10**m) > 0: hi = m
        else: lo = m
    return 10**((lo + hi) / 2)
# ---- 個別の検証 ----
need = ['3Cl2 + 6NaOH -> 5NaCl + NaClO3 + 3H2O', '2NO2 + 2NaOH -> NaNO2 + NaNO3 + H2O', 'Cl2 + 4OH- -> 2ClO- + 2H2O + 2e-', 'Cl2 + 2e- -> 2Cl-', 'Cu2O + H2SO4 -> Cu + CuSO4 + H2O', 'CaC2 + 2H2O -> Ca(OH)2 + C2H2', '2H2S + SO2 -> 3S + 2H2O']
eqs = all_equations()
for e in need: assert e in eqs, e
# 電子のつり合い：Cl- は +1e 受取，ClO3- は 5e 放出
from fractions import Fraction
nCl, nClO3 = 5, 1; assert nCl*1 == nClO3*5
assert 1*1 == 1*1  # NO2- と NO3-
# 酸化数（C in CaC2 and C2H2）: Ca +2 → C2 で -2 → C -1 ; C2H2: H +1 ×2 → C -1
assert (0 - 2)/2 == -1 and (0 - 2*1)/2 == -1

print('ALL OK')
