# 検証スクリプト THE-ION-060  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-060/item.json'
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

lam = dict(H=35, OH=20, Cl=7.5, Na=5.0, Ac=4.0); Ka = 1.8e-5
def cond(hcl, hac, naoh):
    rest = naoh; h_s = max(hcl - rest, 0); rest = max(rest - hcl, 0)
    conv = min(rest, hac); rest -= conv; oh = rest; left = hac - conv
    x = 0
    if oh == 0 and left > 0:
        lo, hi = 0, left
        for _ in range(200):
            x = (lo+hi)/2
            lo, hi = (x, hi) if (h_s + x)*(conv + x) - Ka*(left - x) < 0 else (lo, x)
    return lam['H']*(h_s + x) + lam['Na']*naoh + lam['Cl']*hcl + lam['Ac']*(conv + x) + lam['OH']*oh
c0 = cond(0.010, 0, 0); c1 = cond(0.010, 0, 0.010); c2 = cond(0.010, 0, 0.020)
assert abs(c0 - 0.425) < 1e-12 and abs(c1 - 0.125) < 1e-12 and abs(c2 - 0.375) < 1e-12
assert sig(c1/c0, 2) == 0.29 and sig(c2/c0, 2) == 0.88
sb = (cond(0, 0.010, 0.010) - cond(0, 0.010, 0.009))/0.1   # 1 mL あたり
sa = (cond(0, 0.010, 0.012) - cond(0, 0.010, 0.011))/0.1
print('傾き 前', sb, '後', sa); assert abs(sb - 0.090) < 0.005 and abs(sa - 0.25) < 1e-9
assert sig(0.25/0.090, 2) == 2.8
assert cond(0, 0.010, 0) < 0.05*c0
m = [cond(0.010, 0.010, 0.010*v/10) for v in range(0, 31)]
d1 = [m[i+1]-m[i] for i in range(30)]
assert all(x < 0 for x in d1[:9]) and all(0 < x < 0.012 for x in d1[11:20]) and all(x > 0.02 for x in d1[20:])
print('J156 OK')
