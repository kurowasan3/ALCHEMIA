# 検証スクリプト THE-EQL-059  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-EQL-059/item.json'
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

from sympy import Rational as R, symbols, solve, nsolve
Rg = 8.31e3
pI2 = 0.0100*Rg*360/5.0; pH2 = 0.020*Rg*360/5.0
assert pI2 > 3.0e3
assert sig(pH2 + 3.0e3, 2) == 1.5e4
I2 = 0.0100*0.060; y = 0.0100 - I2
K = (2*y)**2/((0.020-y)*I2); print('K =', K); assert sig(K,2) == 56
# (3) 体積 V で K = n_HI^2/(n_H2 n_I2) となり V によらない → 圧縮後も物質量同じ，濃度 2 倍
yy = symbols('yy')
sols = [s for s in solve(4*yy**2 - 56*(R(1,100)-yy)*(R(4,100)-yy), yy) if 0 < s < R(1,100)]
assert len(sols) == 1
rem = (0.0100 - float(sols[0]))/0.0100; print('残る割合', rem); assert sig(rem*100,2) == 2.3
# 解説中の2次方程式の係数
from sympy import expand
print(expand(4*yy**2 - 56*(R(1,100)-yy)*(R(4,100)-yy)))
print('J132 OK')
