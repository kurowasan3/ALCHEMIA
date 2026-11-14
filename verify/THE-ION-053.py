# 検証スクリプト THE-ION-053  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-053/item.json'
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

from sympy import Rational as R, sqrt, symbols, solve, nsolve
KAgCl, KAg2, KBa = 1.8e-10, 1.2e-12, 1.2e-10
Ag0, Ba0, Cl_B, Cr_B = 2.0e-5, 1.0e-5, 1.0e-3, 5.0e-4
cCl = KAgCl/Ag0; cCr_Ag = KAg2/Ag0**2; cCr_Ba = KBa/Ba0
assert (sig(cCl,2), sig(cCr_Ag,2), sig(cCr_Ba,2)) == (9.0e-6, 3.0e-3, 1.2e-5)
v = lambda c, cB: c*100/cB
vAgCl, vBa, vAg2 = v(cCl, Cl_B), v(cCr_Ba, Cr_B), v(cCr_Ag, Cr_B)
assert (sig(vAgCl,2), sig(vBa,2), sig(vAg2,2)) == (0.90, 2.4, 600)
assert vAgCl < vBa < vAg2
# 0.90 mL で他の塩は未飽和
Cr = Cr_B*vAgCl/100
assert Ba0*Cr < KBa and Ag0**2*Cr < KAg2
# 2.4 mL: AgCl 沈殿後の [Ag+]（体積 100 mL 一定）
Clt = Cl_B*vBa/100
x = symbols('x', positive=True)  # [Ag+]
sol = [s for s in solve(x - KAgCl/x - (Ag0 - Clt), x) if s > 0][0]
Ag = float(sol); Cr = Cr_B*vBa/100
print('2.4 mL での [Ag+] =', Ag, ' Ag2CrO4 のイオン積 =', Ag**2*Cr)
assert Ag**2*Cr < KAg2/100
# 体積変化を考慮した厳密な開始体積でも順序は同じ
import math
def onset(f):
    lo, hi = 1e-6, 1000
    for _ in range(200):
        m = (lo+hi)/2
        lo, hi = (m, hi) if f(m) < 0 else (lo, m)
    return lo
eAgCl = onset(lambda v: Ag0*100/(100+v)*Cl_B*v/(100+v) - KAgCl)
eBa = onset(lambda v: Ba0*100/(100+v)*Cr_B*v/(100+v) - KBa)
print('体積変化込み', eAgCl, eBa); assert eAgCl < eBa
# (5)
v5 = v(cCr_Ba, 1.0e-2); assert sig(v5,2) == 0.12 and v5 < vAgCl
assert sig(v(cCr_Ag, 1.0e-2),2) == 30
print('J131 OK')
