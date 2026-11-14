# 検証スクリプト THE-SOL-053  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-SOL-053/item.json'
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

import sympy as sp, math
D, Vw = sp.Rational(4), 50
r1 = lambda v: 1 / (1 + D*sp.Rational(v)/Vw)
v1 = sp.solve(sp.Eq(1 + D*sp.Symbol('v')/Vw, 100))[0]
assert v1 == sp.Rational(2475, 2) and sig(float(v1), 3) == 1240
e = float(r1(20)**2); u = float(r1(30)*r1(10))
assert sig(e*100, 3) == 14.8 and sig(u*100, 3) == 16.3
best = min(range(0, 41), key=lambda a: float(r1(a)*r1(40 - a)))
assert best == 20
n = 1
while 1.4**n < 100: n += 1
assert n == 14 and n*5 == 70
assert 1.4**13 < 100 <= 1.4**14
assert abs(2/0.146 - 13.7) < 0.05
for V in (40, 60, 80):
    x = 0.08*V; prev = 1
    for k in range(1, 200):
        rn = (1 + x/k)**(-k)
        assert rn < prev + 1e-15 and rn > math.exp(-x)
        prev = rn
Vmin = 12.5*math.log(100)
assert sig(Vmin, 3) == 57.6 and sig(12.5*2*2.303, 3) == 57.6
assert sig(12.5*math.log(1000), 2) == 86 and abs(math.log(1000)/math.log(100) - 1.5) < 1e-12
print('V1', float(v1), '等分', e, '30/10', u, 'n', n, 'Vmin', Vmin)

print('ALL OK')
