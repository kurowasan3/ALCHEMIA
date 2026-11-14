# 検証スクリプト THE-ION-063  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-ION-063/item.json'
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

import math
from scipy.optimize import brentq
K = 10**-6.10; C0, B0 = 1.2e-3, 2.4e-2
def pH_closed(acid):
    # 炭酸の総量 C0+B0 一定，Na+ = B0（陽イオン），Cl- = acid．電荷: h + Na = HCO3- + Cl-
    T = C0 + B0
    f = lambda lh: (lambda h: h + B0 - acid - T*K/(h+K))(10**lh)
    return -brentq(f, -12, -2)
def pH_open(acid):
    # [CO2]=C0 一定: HCO3- = K C0/h，電荷 h + B0 = K C0/h + acid
    f = lambda lh: (lambda h: h + B0 - acid - K*C0/h)(10**lh)
    return -brentq(f, -12, -2)
p0 = pH_closed(0); assert round(p0, 2) == 7.40 and round(pH_open(0), 2) == 7.40
p2 = pH_closed(0.010); p3 = pH_open(0.010); print(p0, p2, p3)
assert round(p2, 2) == 6.20 and round(p3, 2) == 7.17
assert abs(6.10 + (0.70-0.60) - 6.20) < 1e-9 and abs(6.10 + (0.85+0.70-0.48) - 7.17) < 1e-9
xc = brentq(lambda a: pH_closed(a) - 7.10, 0, 0.02); xo = brentq(lambda a: pH_open(a) - 7.10, 0, 0.0239)
print(xc, xo, xo/xc); assert sig(xc,2) == 1.1e-3 and sig(xo,2) == 1.2e-2 and round(xo/xc) == 11
c5 = 0.014/20; assert sig(c5/C0, 2) == 0.58
print('K155 OK')
