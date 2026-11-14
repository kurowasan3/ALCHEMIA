# 検証スクリプト THE-ION-056  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-056/item.json'
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
    keys = ['src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig']
    assert all(k in it for k in keys), [k for k in keys if k not in it]
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R, symbols, solve, sqrt, log, expand
K = R(12,1000)
y = symbols('y')
eq1 = expand((R(2,100)+y)*y - K*(R(2,100)-y)); print(eq1)
assert eq1 == y**2 + R(32,1000)*y - R(24,100000)
y1 = [s for s in solve(eq1, y) if s > 0][0]; h1 = R(2,100)+y1
assert sig(float(y1),2) == 6.3e-3 and sig(float(h1),2) == 2.6e-2 and round(-float(log(h1,10)),1) == 1.6
assert sig(float(y1/R(2,100))*100, 2) == 31
assert abs(float(sqrt(R(1984,10**6))) - 4.454e-2) < 1e-5
z = symbols('z')
eq3 = expand((R(4,100)-z)*(R(12,100)-z) - K*z); assert eq3 == z**2 - R(172,1000)*z + R(48,10000)
z3 = [s for s in solve(eq3, z) if 0 < s < R(4,100)][0]; h3 = R(4,100) - z3
assert sig(float(h3),2) == 5.0e-3 and round(-float(log(h3,10)),1) == 2.3
assert abs(float(sqrt(R(10384,10**6))) - 0.1019) < 1e-4
# (4b) HCl 0.0263 + SO4 0.10
H0, S0 = float(h1), 0.10
w = symbols('w')
wz = [s for s in solve((H0-w)*(S0-w) - float(K)*w, w) if 0 < s < H0][0]
assert H0 - wz < H0; print('(4b) [H+]', H0 - wz)
print('J147 OK')
