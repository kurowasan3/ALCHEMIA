# 検証スクリプト THE-ION-055  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-055/item.json'
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
Ka, Kw = 2.0e-5, 1e-14
def H(Cl, Na, Ac):  # 混合後の濃度
    f = lambda lh: (lambda h: h + Na - Cl - Ac*Ka/(h+Ka) - Kw/h)(10**lh)
    return 10**brentq(f, -14, 0.5)
mix = {'A': (0.05, 0, 0), 'B': (0.05, 0.05, 0.05), 'C': (0.15, 0.05, 0.05), 'D': (0.05, 0.10, 0.10), 'E': (0, 0, 0.05)}
ex = {k: H(*v) for k, v in mix.items()}; print(ex)
ans = {'A': 5.0e-2, 'B': 1.0e-3, 'C': 1.0e-1, 'D': 2.0e-5, 'E': 1.0e-3}
for k in ans: assert abs(ex[k]-ans[k])/ans[k] < 0.015, (k, ex[k])
assert ex['C'] > ex['A'] > ex['B'] and abs(ex['B'] - ex['E'])/ex['E'] < 1e-9 and ex['E'] > ex['D']
print('J145 OK')
