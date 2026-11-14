# 検証スクリプト THE-EQL-058  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-EQL-058/item.json'
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

from sympy import Rational as R, symbols, solve
K1 = R(50,100)/(R(50,100)*1**2); assert K1 == 1
K2 = R(1,4)
Q = lambda h: R(1,2)/(R(1,2)*h**2)
assert Q(1) > K2
q3 = Q(R(3,2)); assert sig(float(q3),2) == 0.44 and q3 > K2
q4 = Q(R(5,2)); assert q4 == R(4,25) and q4 < K2
x = symbols('x', positive=True)
assert solve(Q(1+x) - K2, x) == [1]
assert '0.44' in it['answer'] and '0.16' in it['answer'] and '1.0\\,\\mathrm{mol}' in it['answer']
print('J128 OK')
