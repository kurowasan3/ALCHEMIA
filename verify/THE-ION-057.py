# 検証スクリプト THE-ION-057  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-057/item.json'
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
K1, K2, K3, Kw = 10**-2.1, 10**-7.20, 10**-12.4, 1e-14
def pH(Na, P=0.10):
    def f(lh):
        h = 10**lh; D = h**3 + K1*h**2 + K1*K2*h + K1*K2*K3
        charge_P = P*(K1*h**2 + 2*K1*K2*h + 3*K1*K2*K3)/D
        return h + Na - charge_P - Kw/h
    return -brentq(f, -14, 0)
x = 0.20/3
p2 = pH(0.10 + x); p4 = pH(0.10 + x)   # H3PO4 から: NaOH 0.10+x のみ（NaH2PO4 からの場合と Na の総量が同じ）
print(p2, p4); assert abs(p2 - 7.50) < 0.01 and abs(p4 - 7.50) < 0.01
assert sig(x*40, 2) == 2.7 and sig((0.10+x)*40, 2) == 6.7
r1 = 10**-7.5/K1; r3 = K3/10**-7.5; print(r1, r3); assert r1 < 1e-5 and r3 < 2e-5
print('J151 OK')
