# 検証スクリプト THE-ION-062  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-ION-062/item.json'
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

from scipy.optimize import brentq
import math
K1, K2, Kw = 7.2e-4, 1.8e-5, 1e-14
def exact_pH(C, start='X'):
    # X: HF C + NaAc C ; Y: HAc C + NaF C. 電荷収支: Na+ + H+ = F- + Ac- + OH-（どちらも Na = C, F総 = C, Ac総 = C）
    Na = C
    f = lambda lh: (lambda h: Na + h - C*K1/(h+K1) - C*K2/(h+K2) - Kw/h)(10**lh)
    return -brentq(f, -13, -0.5)
assert K1/K2 == 40
r = 6.32/7.32; assert sig(r,2) == 0.86
h = K2*6.32; assert sig(h,2) == 1.1e-4
pH2 = 4 - 0.057; assert round(pH2,1) == 3.9
pXe = exact_pH(0.10); pYe = exact_pH(0.10, 'Y'); print('厳密 pH(C=0.1)', pXe)
assert round(pXe,1) == 3.9 and abs(pXe - pYe) < 1e-12
assert abs(math.sqrt(K1*K2) - 1.138e-4) < 1e-6
# (4)
C = 1.0e-4
hh = (-K2 + math.sqrt(K2**2 + 4*K2*C))/2; print('近似 h', hh)
assert abs(hh - 3.44e-5) < 0.01e-5 and abs(math.sqrt(75.24) - 8.67) < 0.005
pH4 = 5 - 0.537; assert round(pH4,1) == 4.5
pe = exact_pH(C); print('厳密 pH(C=1e-4)', pe); assert round(pe,1) == 4.5
frac = hh/(hh+K1); print('HF 非解離率', frac); assert 0.04 < frac < 0.05
assert abs(math.log10(1.14) - 0.057) < 0.001 and abs(math.log10(3.44) - 0.537) < 0.001
# (2) の近似の前提: h / [Ac-] が小さい
assert 1.14e-4/(0.137*0.10) < 0.01
print('K134 OK')
