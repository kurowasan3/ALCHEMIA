# 検証スクリプト ORG-ARO-045  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-045/item.json'
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
a, h = sp.symbols('a h', positive=True)
n1 = sp.Rational(100, 1000) * sp.Rational(180, 10) / 1000
n2 = sp.Rational(100, 1000) * (sp.Rational(500, 10) - sp.Rational(184, 10)) / 1000
sol = sp.solve([a + 2*h - n1, 2*a + 2*h - n2], [a, h])
A, H = float(sol[a]), float(sol[h])
M = 9*12 + 8*1.0 + 4*16
assert M == 180
assert abs(A - 1.36e-3) < 1e-9 and abs(H - 2.2e-4) < 1e-9
orig = (A + H) * M / 0.400 * 100
assert sig(orig, 3) == 71.1, orig
assert sig(A * M, 3) == 0.245, A * M
frac = H / (A + H) * 100
assert sig(frac, 2) == 14, frac
print('n1', float(n1), 'n2', float(n2), 'a', A, 'h', H, 'orig%', orig, 'frac%', frac)
assert check_equation('C6H4(OCOCH3)COOH + 3NaOH -> C6H4(ONa)COONa + CH3COONa + 2H2O') is None

print('ALL OK')
