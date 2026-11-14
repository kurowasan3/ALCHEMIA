# 検証スクリプト THE-ION-054  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-054/item.json'
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

from sympy import symbols, nsolve
KA, KB, CA, CB, Kw = 1.8e-5, 1.8e-4, 0.20, 0.20, 1e-14
h = symbols('h', positive=True)
exact = float(nsolve(h - KA*CA/(h+KA) - KB*CB/(h+KB) - Kw/h, h, 6e-3))
approx = (KA*CA + KB*CB)**0.5
print('厳密', exact, '近似', approx)
assert abs(exact-approx)/approx < 0.02 and sig(exact,2) == 6.2e-3 and sig(approx,2) == 6.3e-3
assert abs(39.6**0.5 - 6.29) < 0.005 and abs(0.9**0.5 - 0.949) < 0.001
a_mix = KA/6.29e-3; a_alone = (KA/CA)**0.5
assert sig(a_mix,2) == 2.9e-3 and sig(a_alone,2) == 9.5e-3
hB = (KB*CB)**0.5; assert abs(hB - 6.0e-3) < 1e-12
assert sig(6.29e-3/hB, 2) == 1.0
s_ = hB + (KA*CA)**0.5; print('和', s_); assert sig(s_,2) == 7.9e-3 and sig((KA*CA)**0.5,2) == 1.9e-3
print('J134 OK')
